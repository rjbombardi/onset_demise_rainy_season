#!/usr/bin/python
#========================================================================
#  Subroutine that calculates the end date of the rainy season
# for a time series of precipitation
# nyrs   --> integer for the number of years in the input dataset
# tot    --> total number of points for one year of data (365)
# mtot   --> total number of points in the whole precipitation dataset
# jday   --> an array of Julian days
# day    --> an array of days
# month  --> an array of months
# year   --> an array of years
# jstart --> Julian day of the climatological date when the calculation
#            should start
# precip --> a time series of precipitation anomalies (against mean annual daily)
# npass  --> integer for the number of "passes" for the smoothing of the
#            time series of accumulated precipitation anomalies
#========================================================================
import sys
import math
import numpy
#------------------------------------------------------------------------
# Function that caclulates juilan days
#------------------------------------------------------------------------                 ! 
def rainyseason_demise(nyrs,ytot,jday,day,month,year,jstart,precip,sjday,sday,smonth,syear,curve):
    sseries=numpy.zeros((int(ytot/2)))
    mtot=len(precip)
    tmpday=numpy.zeros((mtot))
    tmpjday=numpy.zeros((mtot))
    tmpmonth=numpy.zeros((mtot))
    tmpyear=numpy.zeros((mtot))
    tmpprec=numpy.zeros((mtot))
    tmpjday[:]=jday[::-1]
    tmpday[:]=day[::-1]
    tmpmonth[:]=month[::-1]
    tmpyear[:]=year[::-1]
    tmpprec[:]=precip[::-1]
    yt=-1
    id=numpy.where(tmpjday[:] == jstart)
    for tt in id[0]:
#    for tt in range(0,mtot-5): # -5 to avoid calcualtion with short time series for last year
#------------------------------------------------------------------------
# Starting the calculation of accumulated anomalies in the rainy season
#------------------------------------------------------------------------                 ! 
#        if tmpjday[tt] == jstart:
        if tt < (mtot-5):         # -5 to avoid calcualtion with short time series for last year
           yt=yt+1
           beg=tt
           ned=beg+int(ytot/2)
           if ned <= mtot-1:  # it is not the last year
              ned2=int(ytot/2)
           if ned > mtot-1:
              ned=mtot-1
              ned2=ned-beg
           sseries[:]=0.
           sseries[0:ned2]=numpy.cumsum(tmpprec[beg:beg+ned2])
           curve[yt,:]=sseries[:]
#-------------------------------------------------------------------------
# Calculating onset and demise of the rainy season
#-------------------------------------------------------------------------
           beg=0
           ons=numpy.where(sseries[0:ned2] == sseries[0:ned2].min())
           if len(ons[0]) > 0:
              beg=ons[0][0]+tt+1
           if beg > 0 and beg < ned: 
              sjday[yt]=tmpjday[beg]
              sday[yt]=tmpday[beg]
              smonth[yt]=tmpmonth[beg]
              syear[yt]=tmpyear[beg]
    sday[:]=sday[::-1]
    sjday[:]=sjday[::-1]
    smonth[:]=smonth[::-1]
    syear[:]=syear[::-1]
    return sjday, sday, smonth, syear,curve
#========================================================================
#                             End of subroutine
#========================================================================

