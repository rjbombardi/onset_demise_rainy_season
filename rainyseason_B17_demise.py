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
import numpy
def rainyseason_B17_demise(nyrs,ytot,jday,day,month,year,jstart,precip,npass,sjday,sday,smonth,syear):
    sseries=numpy.zeros((ytot))
    ssmooth=numpy.zeros((ytot))
    dsdt=numpy.zeros((ytot))
    temp=numpy.zeros((ytot))
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
#------------------------------------------------------------------------
# Starting the calculation of accumulated anomalies in the rainy season
#------------------------------------------------------------------------                 ! 
        if tt < mtot-5:         # -5 to avoid calcualtion with short time series for last year
           yt=yt+1
           beg=tt
           ned=beg+ytot
           if ned <= mtot-1: # it is not the last year
              ned2=ytot
           if ned > mtot-1:
              ned=mtot-1
              ned2=ned-beg
           sseries[:]=0.
           ssmooth[:]=0.
           dsdt[:]=0.
           temp[:]=0.
           sseries[0:ned2]=numpy.cumsum(tmpprec[beg:beg+ned2])
           if ned == mtot-1:  #NEW
              fin=len(sseries[0:ned2])
              tmp=numpy.zeros((fin))
              tmp[:]=sseries[0:ned2]
              if 2*fin > ytot:  #If the series just need a small increment
                 sseries[ned2:ytot]=tmp[::-1][0:ytot-ned2]
                 ned2=ytot
              if 2*fin <= ytot: #If the series is smaller than 1/4 year 
                 sseries[ned2:ned2+fin]=tmp[::-1][0:fin]
                 ned2=ned2+fin
#========================================================================
#          Starting Bombardi and Carvalho (2009) adaptation
#========================================================================
#------------------------------------------------------------------------
# Smoothing the time series of accumulated anomalies
#------------------------------------------------------------------------
           ssmooth[:]=sseries[:]
           temp[:]=sseries[:]
           for nt in range(0,npass):
               temp[0]=0.5*(ssmooth[0]+ssmooth[1])
               temp[ned2-1]=0.5*(ssmooth[ned2-2]+ssmooth[ned2-1])
               temp[1:ned2-1]=0.25*ssmooth[0:ned2-2]+0.50*ssmooth[1:ned2-1]+0.25*ssmooth[2:ned2]
               ssmooth[:]=temp[:]
#------------------------------------------------------------------------
# Calculating the first derivative of sseries
#------------------------------------------------------------------------
           dsdt[:]=0.
           dsdt[0]=ssmooth[1]-ssmooth[0]
           dsdt[ned2-1]=ssmooth[ned2-1]-ssmooth[ned2-2]
# The indices below look weird but that's Python for ya. Last index in range is ignored
           dsdt[1:ned2-1]=0.5*(ssmooth[2:ned2]-ssmooth[0:ned2-2])
#-------------------------------------------------------------------------
# Calculating onset and demise of the rainy season
#-------------------------------------------------------------------------
           beg=0
           st=2
           while st < ned2-2 and beg == 0:    # stops at the firt time the contdition is met
                 if dsdt[st-2] < 0.:
                    if dsdt[st-1] < 0.:
                       if dsdt[st] <= 0.:
                          if dsdt[st+1] > 0.:
                             if dsdt[st+2] > 0.:
                                beg=st+tt+1
                 st=st+1
           if beg > 0: 
              if beg <= mtot-1:    #NEW
                 sjday[yt]=tmpjday[beg]
                 sday[yt]=tmpday[beg]
                 smonth[yt]=tmpmonth[beg]
                 syear[yt]=tmpyear[beg]
    sday[:]=sday[::-1]
    sjday[:]=sjday[::-1]
    smonth[:]=smonth[::-1]
    syear[:]=syear[::-1]
    return sjday, sday, smonth, syear
#========================================================================
#                             End of subroutine
#========================================================================

