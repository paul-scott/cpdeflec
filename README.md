# CPDeflec

Coloured-Pattern Deflectometry is a tool for measuring the slope errors of mirror panels used in solar thermal concentrators.  This tool was developed as part of an engineering honours project.  The motivation was the need for a low cost and simple to use method of characterising the optical properties of parabolic dish mirror panels.  It is currently in a usable state, but requires significant work to set up and calibrate a fully working system (the details of which are not particularly well documented).  The source code is released as a reference for those working on deflectometry.

The current implementation has systematic errors in its output (see paper/thesis below for details).  Further work on the project will likely revolve around identifying and eliminating the source of this error.

## Improvement Ideas

### Code

- Construct unit tests.

- As a priority refactor functions: solveprofile, calcpattvecs, solvesurface.

- Pull out parameters into a configuration file.

### Algorithm

- Implement curve fitting routine for identifying surface positions and norms simultaneously, as an alternative to current extrapolation technique.

### General

- Upload benchmark images for use with testing.

- Upload benchmark colour relation file.

## Thesis

[Solar Concentrating Mirror Panels](doc/thesis.pdf)

## Paper

P. Scott and G. Burgess. [Measurement of mirror panels using coloured-pattern deflectometry](http://stg.anu.edu.au/publications/assets/inproc/scott-solarpaces-2010.pdf). In Proceedings of the 16th SolarPACES Conference, 2010. 
