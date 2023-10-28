# HeliolincRR

HeliolincRR is a HelioLinC variant written in Python.  HelioLinC's initial formulation was published in 2018 ([Holman M., Payne M., Blankley P., et al.](https://iopscience.iop.org/article/10.3847/1538-3881/aad69a)) and extended in 2020 ([Eggl, S., Juric, M., Moeyens, J., et al.](https://ui.adsabs.harvard.edu/abs/2020DPS....5221101E/abstract)) and 2022 ([Heinze, A., Eggl, S., Juric, M., et. al](https://ui.adsabs.harvard.edu/abs/2022DPS....5450404H/abstract)).  In early 2022 I implemented a [version of HelioLinC](https://www.benengebreth.org/dynamic-sky/heliolinc-a-variation-in-6d/) that I've been refining ever since.  This code is product of that work.

All of the above implementations of HelioLinC share three basic traits: 
1. Projection of observer-relative observations to asserted heliocentric positions
2. Propagation of n=2 sized tracklets to a reference epoch (or epochs)
3. Clustering of propagated tracklets at the reference epoch(s)

The 2018 Holman implementation clusters in a 2D angular phase space while the 2020 Eggl and 2022 Heinze implementations cluster position and velocity in 6D cartesian phase space.  All three implementations cluster at a single reference epoch.

## So what's different about HeliolincRR?

There are three noteworthy conceptual differences in HeliolincRR:
1. HeliolincRR uses ***two*** reference epochs and the two propagated tracklet position vectors (thus Heliolinc**RR**) at those epochs as the clustering phase space.  This both fully specifies the orbit of an object and avoids the mixed units and scales problem from using position and velocity in phase space. ([more information here](https://www.benengebreth.org/dynamic-sky/heliolinc-rr/))
2. HeliolincRR attempts to find clusters centered around each propagated tracklet (which allows for overlapping clusters and fewer hypothesis tests) rather than finding mutually exclusive clusters in the phase space. ([here's a visual explanation]())
3. HeliolincRR uses a [fast Lambert solver](https://arxiv.org/abs/1403.2705) for orbit estimation from n=2 sized tracklets. ([how does that work?]())

The first two have substantially improved object recovery in my testing.  The 3rd is mostly about convenience and ease of implementation.  **Run on a two week subset of DP0.3 data, HeliolincRR recovers 99.26% of MBAs and 99.79% of TNOs as pure linkages.**

## Installation

While HeliolincRR is a generalized implementation for any observer at any location, this initial release was designed to work "out of the box" with the [Vera C. Rubin](https://rubinobservatory.org/) telescope's [DP0.3 simulation data](https://dp0-3.lsst.io/index.html) on the Rubin Science Platform (RSP).  As such, the example code in this repository currently requires access to the RSP to run.  With that said, HeliolincRR can be installed like so:

```console
pip install git+https://github.com/bengebre/heliolincrr
```

## Example code

The rendered notebook ```HeliolincRR-TNO.ipynb``` in the [examples/]() directory finds TNOs in a two week subset of DP0.3 data.

## Documentation

Documentation for HeliolincRR is available [here]().

## Other HelioLinC implementations

[HelioLinC3D](https://github.com/lsst-dm/heliolinc2), developed first by Siegfried Eggl and now with Ari Heinze, is the official Rubin implementation of HelioLinC written in C++ and is already [finding PHAs](https://www.nytimes.com/2023/08/05/science/space-asteroids-rubin-heliolinc3d.html) in survey data.

## Acknowlegements

Developer and maintainer:
- [Ben Engebreth](https://benengebreth.org/)

Contributors and collaborators:
- [Siegfried Eggl](https://aerospace.illinois.edu/directory/profile/eggl)
- [Ari Heinze](https://astro.washington.edu/people/aren-heinze)

HeliolincRR was developed by Ben Engebreth, but would not have been possible without the feedback and support of Siegfried Eggl and Ari Heinze.  Siegfried and Ari did not, however, contribute any errors that may be present.  All errors in this work belong to Ben Engebreth alone.

## Footnotes
