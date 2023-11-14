import numpy as np
import pandas as pd
import pickle
import shutil
import os
import logging
import multiprocessing as mp
from numpy import linalg as LA
from itertools import chain
from os.path import exists
from astropy.time import Time
from astropy import units as u
from poliastro.iod import izzo
from poliastro.bodies import Sun
from poliastro.twobody import Orbit
from sklearn.neighbors import KDTree
from tqdm import tqdm

#allow both new and old poliastro interfaces
try:
    from poliastro.twobody.propagation.farnocchia import FarnocchiaPropagator
    POLINEW=True
except:
    from poliastro.twobody.propagation import propagate
    POLINEW=False

class HeliolincRR:
    """
    HeliolincRR: Solar system object linking using heliocentric projection and postition vector clustering at two reference epochs
    """

    def __init__(self,radecs,jdates,obs_locs):
        """
        Instantiates the HeliolincRR class with observational data and sets up instance variables.

        Parameters
        ----------
        radecs : array
            Righ ascension and declination (respectively) of observations as an Nx2 array.
        jdates : array
            Julian dates (MJDs + 2400000.5).
        obs_locs : array
            An Nx3 array of heliocentric observer locations at the corresponding observation times.

        Returns
        -------
        None
        """

        #test that observations are sorted ascending
        if not np.all(np.diff(jdates)[1:] >= 0):
            raise Exception("Observations must be sorted chronologically in ascending order")

        #instance variables
        self.bbf = None 
        self.r_ob_unit = None
        self.combs = None
        self.hypos = None
        self.pfiles = None
        self.hypo_links = None

        #use pandas to unique observation times and observer positions
        df = pd.DataFrame(obs_locs,columns=['sun_obs_x','sun_obs_y','sun_obs_z'])
        df.insert(0,'times',jdates)
        df.drop_duplicates(inplace=True)

        #unique observation times [was 'times']
        self.utimes = np.array(df['times'])

        #sun->observer vectors
        self.r_so = np.array(df[['sun_obs_x','sun_obs_y','sun_obs_z']])

        #dict for frame numnber lookup
        utimes_fnums = dict(zip(self.utimes,np.arange(0,len(self.utimes))))

        #create blobs by frame
        bbf_df = pd.DataFrame(radecs,columns=['ra','dec'])
        bbf_df['bbf_ids'] = np.arange(0,len(jdates))
        bbf_df['fnums'] = [utimes_fnums[t] for t in jdates]
        self.bbf = np.array(bbf_df)

        #unit vectors pointing from observer to observation [equitorial]
        self.r_ob_unit = self.eq2cart(self.bbf[:,0],self.bbf[:,1],1)[0]

    #equitorial spherical (ra/dec) to both equitorial cartesian AND ecliptic cartesian
    @staticmethod
    def eq2cart(ra,dec,d=1):
        """
        Convert equitorial spherical coordinates of right ascension and declination to both 
        equitorial and ecliptic cartesian coordinates with a user-defined magnitude.

        Parameters
        ----------
        ra : float or array-like 
            Right ascension in degrees.
        dec : float or array-like
            Declination in degrees.
        d : float, optional
            Input length to scale the output vector to, by default 1 for unit length.


        Returns
        -------
        eq_cart : array-like
            Equatorial cartesian vectors of the specified magnitude.
        ec_cart : array-like
            Ecliptic cartesian vectors of the specified magnitude.
        """

        ra_rad = np.deg2rad(ra)
        dec_rad = np.deg2rad(dec)

        #obliquity of ecliptic rotation: https://ssd.jpl.nasa.gov/astro_par.html
        ooe = 0.40909262968940374 #rad
        rot_ooe = np.matrix([[1,0,0],[0,np.cos(ooe),np.sin(ooe)],[0,-np.sin(ooe),np.cos(ooe)]])

        x = d*np.cos(ra_rad)*np.cos(dec_rad)
        y = d*np.sin(ra_rad)*np.cos(dec_rad)
        z = d*np.sin(dec_rad)

        eq_cart = np.column_stack((x,y,z)) #equitorial
        ec_cart = np.asarray((rot_ooe*np.matrix(eq_cart).T).T) #ecliptic (p. 269 Vallado)

        return eq_cart, ec_cart

    def create_tracklets(self,dts,dpds):
        """
        Generate tracklets for observations within a specified time range of one another and moving at a specified observer sky rate range.

        Parameters
        ----------
        dts : list or tuple (float, float)
            Minimum and maximum time intervals in days for tracklet formation.
            
        dpds : list or tuple (float, float)
            Minimum and maximum observer relative rates of motion on the sky for tracklet formation.

        Returns
        -------
        combs : list of lists
            List of pairs of observations that meet the conditions specified by the input parameters.

        """

        #dpds in radians
        dpds_rad = np.deg2rad(dpds)

        #tracklet store
        combs = []

        #unit vectors pointing from observer to observation [equitorial]
        #self.r_ob_unit = self.eq2cart(self.bbf[:,0],self.bbf[:,1],1)[0]

        #need to implement magnitude checks still
        #vmags = []
        #use_mags = len(vmags)>0

        #precalc rounded mag difference booleans if vmags and ranges supplied [for mag checks]
        #if use_mags:
        #    mag_bools = {}
        #    rmags = np.round(mags,1)
        #    urmags = np.unique(rmags)
        #    for urmag in urmags:
        #        magrange = magpct*urmag
        #        mag_bools[urmag] = (mags >= (urmag-magrange)) & (mags <= (urmag+magrange))

        #inits
        bbf_ids = self.bbf[:,2].astype(int)
        fnums = self.bbf[:,3].astype(int)
        last_frame = -1
        bbf_len = len(bbf_ids)
        seps_bool_init = np.zeros(bbf_len,dtype=bool)
        ftimes = np.array([self.utimes[f] for f in fnums]) #all times


        with np.errstate(invalid='raise'):
            for bbf_id in tqdm(bbf_ids,desc="generating tracklets",smoothing=0):

                #change in frame numbers/times
                if last_frame!=ftimes[bbf_id]:
                    last_frame = ftimes[bbf_id]
                    diffs = np.abs(ftimes[bbf_id] - ftimes)
                    diff_bool = (diffs>=dts[0]) & (diffs<=dts[1])
                    non_frame_bool = (self.bbf[:,3]>self.bbf[bbf_id,3])
                    combo_bool = (diff_bool) & (non_frame_bool)

                    #if use_mags: [for mag checks]
                    #    combo_bool = combo_bool & mag_bools[rmags[bbf_id]]

                    any_combos = np.any(combo_bool)
                    combo_idx = np.where(combo_bool)[0]

                if ~any_combos:
                    continue

                ##seps_bool = copy.copy(seps_bool_init)

                #angle separations between observation unit vectors
                #seps_calc = vg.angle(self.r_ob_unit[combo_idx,:],np.tile(self.r_ob_unit[bbf_id,:],(len(combo_idx),1))) #thought this would be faster; nope

                #seps_calc = np.rad2deg(np.arccos(np.clip(np.dot(self.r_ob_unit[combo_idx,:],self.r_ob_unit[bbf_id,:]),-1,1))) #deg
                #seps_calc = np.arccos(np.clip(np.dot(self.r_ob_unit[combo_idx,:],self.r_ob_unit[bbf_id,:]),-1,1)) #rad (faster)

                #try except to avoid np.clip() [fastest yet!]
                dot = np.dot(self.r_ob_unit[combo_idx,:],self.r_ob_unit[bbf_id,:])
                try:
                    seps_calc = np.arccos(dot)
                except FloatingPointError as e:
                    print(e)
                    seps_calc = np.arccos(np.clip(dot,-1,1))

                seps_calc = seps_calc/diffs[combo_idx] #sep per time

                #seps_calc_bool = (seps_calc>=dpds[0]) & (seps_calc<=dpds[1]) #deg
                seps_calc_bool = (seps_calc>=dpds_rad[0]) & (seps_calc<=dpds_rad[1]) #rad (faster)

                ##seps_bool[combo_idx] = seps_calc_bool
                ##sub_combs = np.where( (combo_bool) & (seps_bool) )[0]

                #sub_combs = combo_idx[np.where(seps_calc_bool)[0]]
                sub_combs = combo_idx[seps_calc_bool]

                for i in sub_combs:
                    combs.append([bbf_id,i])

        self.combs = combs

    #read with pickle
    @staticmethod
    def pload(fname):
        """
        Loads arbitrary data from file.

        Parameters
        ----------
        fname : str
            Filename to be loaded.

        Returns
        -------
        v : arbitrary
            Returns data stored to specified file. 
        """

        f = open(fname,'rb')
        v = pickle.load(f)
        f.close()

        return v

    #save with pickle
    @staticmethod
    def psave(var,fname):
        """
        Saves arbitrary data to file.

        Parameters
        ----------
        var : arbitrary
            The data that will be saved.
        fname : str
            The filename to write the data to.

        Returns
        -------
        None
        """

        f = open(fname,'wb')
        pickle.dump(var, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    @staticmethod
    def geod2heliod(r_so,r_ob,d):
        """
        Return the multiplier for the observer relative unit vector that scales it to 
        a helioentric vector of length 'd'.

        Parameters
        ----------
        r_so : array
            The sun->observer vector at the observation time.
        r_ob : array
            The observer->body vector at the observation time.
        d : float
            The length of the heliocentric vector you want to create.


        Returns
        -------
        mult : float
            The multiplier that scales the observer->body unit vector to the specified heliocentric distance.

        """

        a = np.sum(r_ob**2)
        b = 2*np.sum(r_ob*r_so)
        c = -1*(d**2 - np.sum(r_so**2))
        r = np.roots([a,b,c])

        mult = r[r>0][0]

        return mult

    @staticmethod
    def get_logger(logger_name):
        """
        Get or initialize basic logger for HeliolincRR runs.

        Parameters
        ----------
        logger_name : str
            Name of the logger and file logs are written to.

        Returns
        -------
        logger : logging.Logger
        """

        logger = logging.getLogger(logger_name)

        if not logger.handlers:
            logger.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s',"%Y-%m-%d %H:%M:%S")
            fh = logging.FileHandler(logger_name + '.log')
            logger.addHandler(fh)
            fh.setFormatter(formatter)
            logger.propagate = 0

        return logger

    def propagate(self,hypos,epochs,cores,max_v=1000,min_init_earth_au=0.1,log=True):
        """
        Multithreading wrapper to Lambert solve (initial orbit determination) all tracklets for all user-defined range hypothesis 
        and then propagate to two user-defined reference epochs.  Calls propagator().

        Parameters
        ----------
        hypos : list of tuples (float, float, float)
            List of all range, range-rate and range-rate-rate hypothesis to test with units of AU, AU/day and AU/day^2 respectively.

        epochs : list or tuple (float, float)
            Two reference epochs (times) that all tracklets will be propagated to for eventual clustering.
            These two times should be Julian dates (not MJDs)

        cores : int
            The number of CPU cores to use to propagate.

        max_v : float, optional
            Do no propagate tracklets that Lambert solve with velocities greater than this (km/s).

        min_init_earth_au : float, optional
            Do not store propagated tracklets that have initial heliocentric projected positions from the observer closer than this (AU).
            Initial projections that are too close to the observer can be a result of geometrically improbable observations and result in mixed linkages.

        log : bool, optional
            Write a log file (propagator.log) to report results of propagation (e.g. failed Lambert IOD for a tracklet).

        Returns
        -------
        None
        """

        #convert epochs to time objects
        oepochs = Time(epochs,format='jd')

        #keep hypos around for reference
        self.hypos = hypos

        #remove old log file if it exists
        if exists('propagator.log'):
            os.remove('propagator.log')

        p = mp.Pool(cores) #mp.cpu_count()
        args = (oepochs,max_v,min_init_earth_au,log)
        pooled_results = [p.apply_async(self.propagator, args=(h,) + args) for h in hypos]
        results = [r.get() for r in tqdm(pooled_results, position=0, leave=True, smoothing=0, desc="propagating")]
        p.close()

        self.pfiles = [r for r in results]

    def propagator(self,hypo,epochs,max_v,min_init_earth_au,log):
        """
        Lambert solve (initial orbit determination) all tracklets using a user-defined range hypothesis 
        and then propagate to two user-defined reference epochs.

        Parameters
        ----------
        hypo : tuple (float, float, float)
            Range, range-rate and range-rate-rate hypothesis to test with units of AU, AU/day and AU/day^2 respectively.

        epochs : list or tuple (float, float)
            Two reference epochs (times) that all tracklets will be propagated to for eventual clustering.
            These two times should be Julian dates (not MJDs)

        cores : int
            The number of CPU cores to use to propagate.

        max_v : float, optional
            Do no propagate tracklets that Lambert solve with velocities greater than this (km/s).

        min_init_earth_au : float, optional
            Do not store propagated tracklets that have initial heliocentric projected positions from the observer closer than this (AU).
            Initial projections that are too close to the observer can be a result of geometrically improbable observations and result in mixed linkages.

        Returns
        -------
        f : str
            Filename of propagated tracklets file written for the specified hypothesis and epochs.
        """

        r,rdot,rdotdot = hypo

        #if file exists we already calculated this
        file_params = ('c'+str(len(self.combs)),) + tuple(np.round(np.array((r,rdot,rdotdot,max_v,min_init_earth_au),dtype=float),8))
        f = "_".join(map(str, file_params)) + '.pkl'
        if exists(f):
            return f

        #count times propagated tracklet exceeds max_v 
        max_v_count = 0

        #log
        if log:
            logger = self.get_logger('propagator')

        #times as objects for integrator
        times = [Time(t,format='jd') for t in self.utimes]

        #dts from start
        dts = np.asarray([(t - times[0]).value for t in times])

        #frame numbers (time index) for each bbf
        fnums = self.bbf[:,3].astype(int)

        #results
        prvs = [] #propagated vectors
        porbs = [] #propagated orbits
        pcids = [] #propagated combination ids

        num_tracklets = len(self.combs)
        if log:
            msg = "Starting propagation for hypothesis: ({:.2f},{:.5f},{:.8f}) with {:d} tracklets".format(r,rdot,rdotdot,num_tracklets)
            logger.info(msg)

        for i,c in enumerate(self.combs):
            #combination times
            t0 = times[fnums[c[0]]]
            t1 = times[fnums[c[1]]]

            #earth->object vectors and sun->earth vectors
            r_ob0 = self.r_ob_unit[c[0],:]
            r_ob1 = self.r_ob_unit[c[1],:]
            r_so0 = self.r_so[fnums[c[0]],:]
            r_so1 = self.r_so[fnums[c[1]],:]

            #computed range
            range_dt0 = r + rdot*dts[fnums[c[0]]] + 0.5*rdotdot*dts[fnums[c[0]]]**2
            range_dt1 = r + rdot*dts[fnums[c[1]]] + 0.5*rdotdot*dts[fnums[c[1]]]**2

            #solve for specified helio distance
            try:
                alpha0 = self.geod2heliod(r_so0,r_ob0,range_dt0)
                alpha1 = self.geod2heliod(r_so1,r_ob1,range_dt1)
            except Exception as e:
                if log:
                    msg = "No physical heliocentric projection possible from geometry supplied - H:({:.2f},{:.5f},{:.8f}) T:{:d}/{:d}".format(r,rdot,rdotdot,i,num_tracklets)
                    logger.debug(msg)
                continue

            #heliocentric vector too close to earth [speed up by not calculating here and below]
            if  (alpha0 <= min_init_earth_au) or (alpha1 <= min_init_earth_au):
                if log:
                    msg = "Initial heliocentric projection vector too close to observer - H:({:.2f},{:.5f},{:.8f}) T:{:d}/{:d}".format(r,rdot,rdotdot,i,num_tracklets)
                    logger.debug(msg)
                continue

            #convert earth->object vectors and sun->earth vectors to heliocentric vector of length 'd'
            hr0 = (r_so0 + alpha0*r_ob0)*u.AU
            hr1 = (r_so1 + alpha1*r_ob1)*u.AU

            #calc heliocentric velocity vectors to connect observations (km/sec)
            try:
                velocities = izzo.lambert(Sun.k, hr0, hr1, t1-t0, rtol=1e-4)
                if POLINEW:
                    v0, v1 = velocities
                else:
                    (v0, v1), = velocities

            except Exception as e:
                if log:
                    msg = "Lambert solver can't create orbit from tracklet - H:({:.2f},{:.5f},{:.8f}) T:{:d}/{:d}".format(r,rdot,rdotdot,i,num_tracklets)
                    logger.debug(msg)
                continue

            #if velocity [km/s] not crazy propagate back to comparison epoch
            if LA.norm(v1.value) < max_v:
                orb = Orbit.from_vectors(Sun, hr1, v1)

                #method 1 (original)
                #o0 = orb.propagate(epochs[0]-t1,rtol=1e-8) #propagate to epoch0
                #o1 = orb.propagate(epochs[1]-t1,rtol=1e-8) #propagate to epoch1

                #method 2 of propagation: https://github.com/poliastro/poliastro/issues/1008 (faster)
                timedelta_vector = epochs - t1
                if POLINEW:
                    prv = np.reshape(np.transpose(FarnocchiaPropagator.propagate_many(orb,orb._state,timedelta_vector)[0].to(u.AU).value),(1,6))[0].tolist()
                else:
                    prv = np.reshape(np.transpose(propagate(orb, timedelta_vector, rtol=1e-4).xyz.to(u.AU).value),(1,6))[0].tolist()

                #prvs.append(list(np.concatenate((o0.r.to(u.AU).value,o1.r.to(u.AU).value)))) #for method 1
                if ~np.any(np.isnan(prv)):
                    prvs.append(prv) #for method 2
                    pcids.append(c)
                elif log:
                    msg = "NAN in propagated state - H:({:.2f},{:.5f},{:.8f}) T:{:d}/{:d}".format(r,rdot,rdotdot,i,num_tracklets)
                    logger.debug(msg)
            else:
                max_v_count+=1

        #log how many times propagated velocity exceeded max_v
        if log:
            msg = "Propagated velocity greater than max_v {:d} times - H:({:.2f},{:.5f},{:.8f}) T:{:d}/{:d}".format(max_v_count,r,rdot,rdotdot,num_tracklets,num_tracklets)
            logger.debug(msg)


        #return data to file
        prvs = np.array(prvs)
        porbs = np.array(porbs)
        ret = [prvs, porbs, pcids, hypo]
        f_tmp = f + '.tmp'
        self.psave(ret, f_tmp)
        shutil.move(f_tmp, f)

        return f

    @staticmethod
    def min_nights_cluster_filter(links,ccids,bbf,times,n):
        """
        Filter links that don't have enough nights of data.

        Parameters
        ----------
        links : list of lists
            Links that you are filtering.
        ccids : list of lists
            Combination ids (aka tracklet pair ids).
        bbf : array
            Nx4 array of observations.
        times : array
            Unique observation times.
        n : int
            Minimum number of nights a link must have.

        Returns
        -------
        accepted_links : list of lists
            The links that have >= n nights of observations.
        accepted_ccids : list of lists
            The combination pairs of the links.

        """
        accepted_links = []
        accepted_ccids = []
        for i,link in enumerate(links):
            ltimes = times[np.array([int(b) for b in bbf[link,3]])]
            nights = np.sum(ltimes[1:]-ltimes[0:-1] > 0.5) + 1
            if nights>=n:
                accepted_links.append(link)
                accepted_ccids.append(ccids[i])

        return accepted_links, accepted_ccids

    @staticmethod
    def flatten_links(dlinks):
        """
        Turn links at each hypothesis into one list of links.
        Flattens a list of list of lists into one list of lists.

        Parameters
        ----------
        dlinks : list of list of lists
            List of links for all hypotheses.

        Returns
        -------
        links : list of lists
            Flattened list of list of lists.
        """

        links = []
        for ll in dlinks:
            for l in ll:
                links.append(l)

        return links

    @staticmethod
    def ulinks(links):
        """
        Unique links.

        Parameters
        ----------
        links : list of lists
            Candidate links.

        Returns
        -------
        ulinks : list of lists
            Uniqued candidate links.
        """

        try:
            final_links = [list(l) for l in set(tuple(l) for l in links)]
        except Exception:
            links = HeliolincRR.flatten_links(links)
            final_links = [list(l) for l in set(tuple(l) for l in links)]

        return final_links

    def cluster(self,tol,min_len,min_nights,cores,max_tracklets=30):
        """
        Multithreading wrapper for clustering propagated trackets to generate candidate linkages.

        Parameters
        ----------
        tol : float
            The clustering tolerance in AU.
        min_len : int
            The minimum length (number of observations) of a candidate link.
        min_nights : int
            The minimum number of nights of observations in a candidate link.
        cores : int
            The number of CPU cores to use to propagate.
        max_tracklets : int, optional
            The maximum amount of tracklets in a single candidate link.  This prevents link length from
            blowing up in certain scenarios.

        Returns
        -------
        hypo_links : list of lists
            List of candidate links for each range hypothesis tested.
        """

        p = mp.Pool(cores) #mp.cpu_count()
        args = (tol,min_len,min_nights,max_tracklets)
        pooled_results = [p.apply_async(self.clusterer, args=args + (pfile,)) for pfile in self.pfiles]
        hypo_links = [r.get() for r in tqdm(pooled_results, position=0, leave=True, smoothing=0, desc="clustering")]
        p.close()

        self.hypo_links = hypo_links

        return hypo_links

    def clusterer(self,tol,min_len,min_nights,max_tracklets,fname):
        """
        Cluster propagated trackets to generate candidate linkages.

        Parameters
        ----------
        tol : float
            The clustering tolerance in AU.
        min_len : int
            The minimum length (number of observations) of a candidate link.
        min_nights : int
            The minimum number of nights of observations in a candidate link.
        max_tracklets : int, optional
            The maximum amount of tracklets in a single candidate link.  This prevents link length from
            blowing up in certain scenarios.
        fname : str
            The filename of the propagated tracklet file that will be read.

        Returns
        -------
        links : list of lists
            List of candidate links for each range hypothesis tested.
        """

        #fcombs is list if list of frame numbers for the tracklet combinations
        #combs is list of list of bbf ids of the tracklet combinations

        #use data from file (overriding input drvs/combs)
        if fname:
            #fresults = self.pload(fname)
            drvs_clust, orbs, combs, hypo = self.pload(fname) 
            #drvs_clust = fresults[0]
            #combs = fresults[2]
            #hypo = fresults[3]
            #del(fresults)

        #frame nums for observations
        fnums = self.bbf[:,3].astype(int)

        #no vectors to process
        if (np.shape(drvs_clust)[0])==0:
            return []

        #handle case where there aren't many trancklets
        num_tracklets = len(drvs_clust[:,0])
        if num_tracklets < max_tracklets:
            max_tracklets = num_tracklets

        #frame combinations
        fcombs = [[fnums[c[0]],fnums[c[1]]] for c in combs]

        #normalize phase space
        #drvs_clust = copy.deepcopy(drvs) #change to copy.copy?
        #mean_d0 = np.mean(LA.norm(drvs_clust[:,0:3],axis=1))
        #mean_d1 = np.mean(LA.norm(drvs_clust[:,3:6],axis=1))
        drvs_clust[:,0:3] = drvs_clust[:,0:3]/hypo[0] #previously mean_d0
        drvs_clust[:,3:6] = drvs_clust[:,3:6]/hypo[0] #previously mean_d1
        #del(drvs)

        #build KDTree
        tree = KDTree(drvs_clust)

        #init
        links = [] #final list of links
        ccids = [] #final list of clustered bbf_id combinations

        for i,rv in enumerate(drvs_clust):

            dist, kidx = tree.query(np.asarray([rv]),k=max_tracklets)
            dist = dist[0]
            kidx = kidx[0]

            #tolerance test
            kidx = kidx[np.where(dist<=tol)]

            if len(kidx)>1:
                ccs = [combs[i]] #clustered bbf_id combs
                cfcombs = [fcombs[i]] #clustered frame number combs

                #if new tracklet adds the same number of frames as bbf_ids then keep
                #if it doesn't it means there are two different detections in the same frame
                for j,k in enumerate(kidx):
                    new_ccs = ccs + [combs[k]]
                    new_cfcombs = cfcombs + [fcombs[k]]
                    #if len(np.unique(new_ccs)) == len(np.unique(new_cfcombs)): #OLD
                    if len(set(list(chain.from_iterable(new_ccs)))) == len(set(list(chain.from_iterable(new_cfcombs)))): #OPTIMIZED
                        ccs = new_ccs
                        cfcombs = new_cfcombs

                #link = list(np.unique(ccs)) #link ids from combination (bbf) ids
                link = sorted(set(list(chain.from_iterable(ccs)))) #link ids from combination (bbf) ids OPTIMIZED?

                #uniq fnum length long enough and no dupe frames
                if len(link) >= min_len:
                    links.append(link)
                    ccids.append(ccs)

        #filter by min nights
        links = self.min_nights_cluster_filter(links,ccids,self.bbf,self.utimes,min_nights)[0]

        #unqiue
        links = self.ulinks(links)

        return links
