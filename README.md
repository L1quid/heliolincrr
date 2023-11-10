# HeliolincRR

HeliolincRR is a HelioLinC variant for solar system object linking written in Python.  HelioLinC's original formulation was published in 2018 ([Holman M., Payne M., Blankley P., et al.](https://iopscience.iop.org/article/10.3847/1538-3881/aad69a)) and extended in 2020 ([Eggl, S., Juric, M., Moeyens, J., et al.](https://ui.adsabs.harvard.edu/abs/2020DPS....5221101E/abstract)) and 2022 ([Heinze, A., Eggl, S., Juric, M., et. al](https://ui.adsabs.harvard.edu/abs/2022DPS....5450404H/abstract)).  In early 2022 I implemented a [version of HelioLinC](https://www.benengebreth.org/dynamic-sky/heliolinc-a-variation-in-6d/) that I've been refining ever since.  This code is product of that work.

All of the above implementations of HelioLinC share three basic traits: 
1. Projection of observer-relative observations to asserted heliocentric positions
2. Propagation of n=2 sized tracklets to a reference epoch (or epochs)
3. Clustering of propagated tracklets at the reference epoch(s)

The 2018 Holman implementation clusters in a 4D parameterized angular phase space while the 2020 Eggl and 2022 Heinze implementations cluster position and velocity in 6D cartesian phase space.  All three implementations cluster at a single reference epoch.

## So what's different about HeliolincRR?

There are three noteworthy conceptual differences in HeliolincRR compared to the implementations above:
1. HeliolincRR uses ***two*** reference epochs and the two propagated tracklet position vectors at those epochs as the clustering phase space (thus Heliolinc**RR**).  This both fully specifies the orbit of an object and avoids the mixed units and scale problem from using position and velocity in phase space. ([more information here](https://www.benengebreth.org/dynamic-sky/heliolinc-rr/))
2. HeliolincRR attempts to find clusters centered around each propagated tracklet (which allows any single tracklet to belong to multiple clusters per hypothesis) rather than finding only mutually exclusive clusters in each hypothesis tested. ([here's a visual explanation](#f1))
3. HeliolincRR uses a [fast Lambert solver](https://arxiv.org/abs/1403.2705) for orbit determination of n=2 sized tracklets. ([what's a Lambert solver?](#f2))

The first two differences have substantially improved object recovery in my testing compared to *my* prior single epoch position and velocity phase space implementation.  **HeliolincRR recovers 99.26% of MBAs and 99.79% of TNOs as pure linkages from the set of *all* solar system objects in the DP0.3 data for an arbitrary two week interval.**  

The 3rd difference is mostly about convenience and ease of implementation.  A Lambert solver can quickly find the Keplerian orbit (elliptical, parabolic or hyperbolic) that joins two position vectors at two observation times, which is exactly what heliocentric projected tracklets are.  Noteably, the formulation of Lambert's problem as two position vectors at two observation times uniquely specifying an orbit is what led to the insight of using two position vectors at two reference epochs for the clustering phase space.

## Installation

While HeliolincRR is a generalized implementation for any observer at any location, this initial release was designed to work "out of the box" with the [Vera C. Rubin](https://rubinobservatory.org/) telescope's [DP0.3 simulation data](https://dp0-3.lsst.io/index.html) on the Rubin Science Platform (RSP).  As such, the example code in this repository currently requires access to the RSP to run.  With that said, the code should be explanatory for those looking to use HeliolincRR with observations from other observatories.  

HeliolincRR can be installed like so:

```console
pip install git+https://github.com/bengebre/heliolincrr
```

## Example code

The rendered notebook ```HeliolincRR-TNO.ipynb``` in the [examples/Rubin-LSST](https://github.com/bengebre/heliolincrr/tree/main/examples/Rubin-LSST) directory finds TNOs in a two week subset of DP0.3 data.

## Other HelioLinC implementations

[Heliolinc3D](https://github.com/lsst-dm/heliolinc2), developed first by Siegfried Eggl and now Ari Heinze, is the official Rubin implementation of HelioLinC written in C++ and is already [finding PHAs](https://www.nytimes.com/2023/08/05/science/space-asteroids-rubin-heliolinc3d.html) in survey data.

## Acknowlegements

Developer and maintainer:
- [Ben Engebreth](https://benengebreth.org/)

Contributors and collaborators:
- [Siegfried Eggl](https://aerospace.illinois.edu/directory/profile/eggl)
- [Ari Heinze](https://astro.washington.edu/people/aren-heinze)

HeliolincRR was developed by Ben Engebreth, but would not have been possible without the feedback and support of Siegfried Eggl and Ari Heinze.  Siegfried and Ari did not, however, contribute any errors that may be present in HeliolincRR.  All errors in this work belong to Ben Engebreth alone.

HeliolincRR uses [poliastro](https://github.com/poliastro/poliastro) for orbit propagation and a fast (Izzo) Lambert solver implementation.  

[DP0.3](https://dp0-3.lsst.io/index.html) is a Data Preview of simulated solar system objects for the [Rubin Observatory's](https://rubinobservatory.org/) [Legacy Survey of Space and Time (LSST)](https://rubinobservatory.org/explore/lsst) created by the [Solar System Science Collaboration](https://lsst-sssc.github.io/).

Specifically, The DP0.3 data set was generated by members of the Rubin Solar System Pipelines and Commissioning teams, with help from the LSST Solar System Science Collaboration, in particular: Pedro Bernardinelli, Jake Kurlander, Joachim Moeyens, Samuel Cornwall, Ari Heinze, Steph Merritt, Lynne Jones, Siegfried Eggl, Meg Schwamb, Grigori Fedorets, and Mario Juric.


## Footnotes

<a name="f1">1</a>. A mutually exclusive clustering algorithm could reasonably choose the clusters shown in red circles in the left image for this synthetic tracklet phase space.  But both of those clusters contain mixed data (blue in the gold set and gold in the blue set).  Therefore, a pure cluster cannot be recovered.  However, a clustering algorithm that attempts to form a cluster around each tracklet (right image), can find numerous pure clusters with the same clustering radius for the same data.  Green circles indicate pure clusters and red circles indicate impure clusters.

![mutually exclusive clusters](https://benengebreth.org/misc/me.png?2)
![overlapping clusters](https://benengebreth.org/misc/ol.png?2)

<a name="f2">2</a>. Lambert solvers find orbits given two position vectors and a time of flight between those two position vectors.  Because all versions of HelioLinC project angular observations to heliocentric positions, every n=2 sized tracklet is two position vectors at two observation times.  A [fast Lambert solver](https://arxiv.org/abs/1403.2705) can thus determine the orbit that connects the two asserted heliocentric positions which can then be propagated to a reference epoch or epochs.  See also: [Lambert's problem](https://en.wikipedia.org/wiki/Lambert%27s_problem).
