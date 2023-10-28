# HeliolincRR

HeliolincRR is a HelioLinC variant written in Python.  HelioLinC's initial formulation was published in 2018 (Holman M, Payne M, Blankley P., et al.) and extended in 2020 (Eggl, S., Juric, M., Moeyens, J., et al.) and 2022 (Heinze, A., Eggl, S., Juric, M., et. al).  In early 2022 I implemented a version of HelioLinC that I've been refining ever since.  This code is product of that work.

All of the above implementations of HelioLinC share three basic traits: 
1. Projection of observer-relative equatorial coordinate observations to asserted heliocentric positions
2. Propagation of n=2 sized tracklets to a reference epoch (or epochs)
3. Clustering of tracklets at the reference epoch(s)

The 2018 Holman implementation clusters in a 2D angular phase space while the 2020 Eggl and 2022 Heinze implementations cluster position and velocity in 6D cartesian phase space.  All three implementations cluster at a single reference epoch.

## So what's different about HeliolincRR?

There are three noteworthy conceptual differences in HeliolincRR:
1. HeliolincRR uses **two** reference epochs and the two position vectors at those epochs as the clustering phase space ([why?](https://www.benengebreth.org/dynamic-sky/heliolinc-rr/))
2. HeliolincRR attempts to find clusters centered around each propagated tracklet (which allows for overlapping clusters) rather than finding mutually exclusive clusters in the phase space. ([tell me more]()]
3. HeliolincRR uses a [fast Lambert solver]() for orbit estimation from n=2 sized tracklets. ([how's that work?]()]

