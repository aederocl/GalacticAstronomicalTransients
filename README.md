# GalacticAstronomicalTransients

## Description:
This script queries VSX (via VizieR), Simbad, Gaia and WISE to find counterparts
of transients, given the coordinates of an object and a search radius.

## Dependencies:
- numpy
- matplotlib
- astropy
- astroquery
- pyvo

## usage: 
checkGAT.py [-h] -f F -r R ra dec

## Examples:
  python checkGAT.py -f d -r 5 303.8972306 -14.2790126

  python checkGAT.py -f x -r 5 07:25:53.03 +29:04:10.2

## Todo:
- add GALEX
- add eROSITA
- add light curve plot: ASAS-SN
- add J-surveys
- make more elegant to avoid plotting in case no objects are found
