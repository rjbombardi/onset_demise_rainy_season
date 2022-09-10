#!/usr/bin/python
"""
Main program that calculates the characteristics of the rainy and dry seasons
"""
import sys
import numpy as np
from scipy.stats import norm
import netCDF4 as nc
import math
from datetime import timedelta, date
# importing local functions
sys.path.append('/home/rodrigo/python_programs/functions/') #Edit this path
from rainyseason_onset import rainyseason_onset
from rainyseason_B17_onset import rainyseason_B17_onset
from rainyseason_demise import rainyseason_demise
from rainyseason_B17_demise import rainyseason_B17_demise
"""
Program that calculates the characteristics of the rainy and dry seasons:
onset and demise dates; duration; and accumualted precipitation

The mean annual cycle is smoothed using the first three harmonics of the mean annual
cycle

The input variables below need to be defined by the user. Either retrieved from the
dataset or created by the user.

The program saves the output as NetCDF files.

Input:
  yr0    : first year of data
  tot    : total number of data in a single year (can be retrieved from dataset)
  ntot   : total number of data points in time (can be retrieved from dataset)
  nlat   : total number of data points in latitude (can be retrieved from dataset)
  nlon   : total number of data points in longitude (can be retrieved from dataset)
  lats   : array with latitude values (can be retrieved from dataset)
  lons   : array with longitude values (can be retrieved from dataset)
  missval: value for missing values (can be retrieved from dataset)
  dper   : minimum percentage [0.,100.] of missing data that can be tolerated
  pathout: path of the directory to where the output will be saved
  prefix : The first part of the name of output files.
Output:
  A set of NetCDF files containing gridded values
  onset_jday   : onset date in Julian days or Day of Year
  onset_day    : day of the month the onset happened
  onset_month  : month the onset happened
  onset_year   : year the onset happened
  demise_jday  : demise date in Julian days or Day of Year
  demise_day   : day of the month the demise happened
  demise_month : month the demise happened
  demise_year  : year the demise happened
  totwet       : accumulated precipitatin during the rainy season
  totdry       : accumulated precipitatin during the dry season
  durwet       : duration of the rainy season
  durdry       : duration of the dry season

"""
#====================================== Reading Data ===================================
"""
This block needs to be edited by the user. Read your input dataset and retrieve (or create)
the Input variables described above.
"""

dper=25.
yr0=
ntot=
nlat=
nlon=

""" Example of how to create latitude and longitude arrays"""
lats=np.arange(-90,90,0.5)
lons=np.arange(0.,360.,0.5)

""" Example of how to create the dataset array"""
prec=np.zeros((ntot,nlat,nlon))

"""Example of how to find the value for missing data. Usually a very large nagative or positive number. For a very large positive number, change .min() to .max() """
missval=prec.min()   # missing value

#=======================================================================================
#---------------------------------------------------------------------------------------
"""                  No Further Editing is Required from this point on               """
#---------------------------------------------------------------------------------------
#=======================================================================================
print("Data read.")

def julian(dd,mm,yy,noleap=1):
    """
    Function that calculates julian days [Day of Year] from day,month, and year imput
    Imput:
       yy:     year [integer]
       mm:     month [integer]
       dd:     day [integer]
       noleap = flag to indicate whether or not leap years should be considered.
       noleap = 0 --> time searies contain Feb 29
       noleap = 1 --> time series does not contain Feb 29
   Output:
       jday:   COrresponding Julian day or Day of Year [1,366]
    Example:
    --------
      >>> jday = julian(1,1,1980) # no leap years
      >>> jday = julian(1,1,1980,0)
    """
    mon=[31,28,31,30,31,30,31,31,30,31,30,31]
    if noleap==0:
       if yy % 4 == 0 and yy % 100 != 0 or yy % 400 == 0:
          mon=[31,29,31,30,31,30,31,31,30,31,30,31]
    if mm == 1:
       jday=dd
    if mm > 1:
       jday=sum(mon[0:mm-1])+dd
    return jday

#=======================================================================================
def dates_daily(y0,m0,d0,mtot,noleap):
    """
    Function that calculates arrays of dates (hour, day,month, year)

    Imput:
       y0:     initial year [integer]
       m0:     initial month [integer]
       d0:     initial day [integer]
       mtot:   total number of elements in the time series [integer]
       noleap: flag to indicate whether or not leap years should be considered.
               noleap = 0 --> time searies contain Feb 29
               noleap = 1 --> time series does not contain Feb 29
    Output:
       day:    array with values of days
       month:  array with values of months
       year:   array with values of years
    Example:
    --------
      >>> day,month,year=dates_daily(1980,1,1,365,0)
    """
    day=np.zeros((mtot))
    month=np.zeros((mtot))
    year=np.zeros((mtot))
    start_date=date(y0,m0,d0)
    deltad = timedelta(days=1)
    single_date=start_date
    dt=0
    while dt < mtot:
          day[dt]=single_date.day
          month[dt]=single_date.month
          year[dt]=single_date.year
          if noleap == 0:
             single_date=single_date+deltad
          if noleap == 1:
             if year[dt] % 4 != 0:
                single_date=single_date+deltad
             if year[dt] % 4 == 0:
                if year[dt] % 100 != 0:
                   if month[dt] != 2:
                      single_date=single_date+deltad
                   if month[dt] == 2:
                      if day[dt] < 28:
                         single_date=single_date+deltad
                      if day[dt] == 28:
                         single_date=single_date+2*deltad
                if (year[dt] % 100 == 0) & (year[dt] % 400 != 0):
                   single_date=single_date+deltad
                if year[dt] % 400 == 0:
                   if month[dt] != 2:
                      single_date=single_date+deltad
                   if month[dt] == 2:
                      if day[dt] < 28:
                         single_date=single_date+deltad
                      if day[dt] == 28:
                         single_date=single_date+2*deltad
          dt=dt+1
    return day, month, year

#=======================================================================================
"""
Funtion that calculates the Fourier coefficients and the explained variance of the Nth
first harmonics of a time series

Input:
   tseries: input time series
   nmodes : number of harmonics to retain (N)
   coefa  : Array with N (or 'nmodes') elements
   coefb  : Array with N (or 'nmodes') elements
   hvar   : Array with N (or 'nmodes') elements
   missval: Falg value for missing data
Output:
   coefa: Array of A coefficients of the Nth first harmonics
   coefb: Array of B coefficients of the Nth first harmonics
   hvar : Array of explained variance of the Nth first harmonics
"""
def Harmonics(coeafa,coefb,hvar,tseries,nmodes,missval):
    mtot=len(tseries)               # retrieving the lenght of the time dimension
    time=np.arange(1,mtot+1,1.)     # Just a an array of increasing numbers
    tdata=tseries[:]                # Adjusting the mean annual cycle
    tdata[tseries==missval]=0.      # Padding missing values with zeros just to be safe
    svar=sum((tdata[:]-np.mean(tdata))**2)/(mtot-1)
    nm=nmodes
    if 2*nm > newdim:
       nm=mtot/2
    coefa=np.zeros((nm))
    coefb=np.zeros((nm))
    hvar=np.zeros((nm))
    for tt in range(0,nm):
        Ak=np.sum(tdata[:]*np.cos(2.*math.pi*(tt+1)*time[:]/float(mtot)))
        Bk=np.sum(tdata[:]*np.sin(2.*math.pi*(tt+1)*time[:]/float(mtot)))
        coefa[tt]=Ak*2./float(mtot)
        coefb[tt]=Bk*2./float(mtot)
        hvar[tt]=mtot*(coefa[tt]**2+coefb[tt]**2)/(2.*(mtot-1)*svar)
    return coefa,coefb,hvar
  
#================================= Formatting Data ======================================

#------------------------------------------------------------------------
# Creating vectors of dates
#------------------------------------------------------------------------                 !
yr0=1979
day,month,year=dates_daily(yr0,1,1,ntot,0)

print('data read')
print('Formating data...')

#------------------------------------------------------------------------
# FInding missing value
#------------------------------------------------------------------------
missval=-999.000   # missing value
#================================= Formatting Data ======================================
print("Formatting Data...")

"""
Removing Feb 29th. This block averages Feb 28 and 29 in leap years. 
"""
id=np.where((month == 2.) & (day == 29.))
id2=[x-1 for x in id]  # creating a list equal to (id-1) and python is stupid
prec[id2[0],:,:]=0.5*(prec[id[0],:,:]+prec[id2[0],:,:])
prec=np.delete(prec,id,axis=0)
year=np.delete(year,id,axis=0)
month=np.delete(month,id,axis=0)
day=np.delete(day,id,axis=0)
ntot=len(year)
prec[prec<0.]=missval

#------------------------------------------------------------------------
# Function that caclulates juilan days
#------------------------------------------------------------------------                 !
jday=np.zeros((ntot))
for tt in range(0,ntot):
    jday[tt]=julian(int(day[tt]),int(month[tt]),int(year[tt]))

#---------------------------------------------------------------------------------------
"""
Masking Missing values. Sometimes datasets have significant amounts of
missing data (e.g. land only data, ocean only data, complex topography).
In such cases it is sometimes useful to mask those regions. Masking can
prevent code errors and improve code efficiency. The minimum percentage
of missing data allowed is determined by namelist variable "dper". If a
grid pint has more missing data then the minimum percentage, that grid
point will be masked at all times.
"""
thres=(ntot-0.01*dper*ntot) # minimum threshold of non-missing data
mask=np.zeros((nlat,nlon))
for it in range(0,nlat):
    for jt in range(0,nlon):
        id=np.where(prec[:,it,jt]!=missval)
        if len(id[0]) >= thres:
           mask[it,jt]=1.

print("Data Formatted.")

#=======================================================================================
"""

Computing input information for the calculation of the characteristics of the rainy
and dry seasons

"""
#=======================================================================================
"""
Calculating the daily annual mean [mm/day]

"""
rm=np.zeros((nlat,nlon))
for it in range(0,nlat):
    for jt in range(0,nlon):
        if mask[it,jt] == 1.:  # removing grid points that have more than 25% of missing data
           tmp=prec[:,it,jt]
           id=np.where(tmp >= 0.)  # removing missing data
           if len(id[0]) > 1:      # have to specify id[0] because id is a tuple
              rm[it,jt]=np.mean(tmp[id[0]])

print("Calculating the mean annual cycle... This might take a while.")
"""
This block will calculate the mean annual cycle for the whole time series. It will need
to be adapted if the period of interest is only a portion of the total time series.
For example, climatologies of a 30 year period of reference: 1981-2010.
"""

""" calculating the mean annual cycle for each grid point """
cycle=np.zeros((tot,nlat,nlon))
for tt in range(0,tot):
    id=np.where(jday[:] == jday[tt])
    for it in range(0,nlat):
        for jt in range(0,nlon):
            if mask[it,jt] == 1.:
               tmp=prec[id[0],it,jt]
               id2=np.where(tmp >= 0.)
               if len(id2[0]) > 1:      # have to specify id2[0] because id2 is a tuple
                  cycle[tt,it,jt]=np.mean(tmp[id2[0]])


print("Mean annual cycle calculated.")

print("Smoothing the mean annual cycle... This might take a while.")
"""
This block will calculate the smoothed mean annual cycle. This is anessential part of
calculating anomalies (often overlooked).

For details, see: Bombardi RJ and Carvalho LMV (2017). Simple Practices in Climatological Analyses: A Review. Revista Brasileira de Meteorologia. 32 (3), 311-320

This block also calculates the explained variance of the first 3 harmonics of the
mean annual cycle. This is used to mask regions with zero or multiple rainy seasons

"""
time=None
time=np.arange(1,tot+1,1.)
harm1=np.zeros((nlat,nlon))
harmonic1=np.zeros((tot,nlat,nlon))
harm2=np.zeros((nlat,nlon))
harm3=np.zeros((nlat,nlon))
smoothed=np.zeros((tot,nlat,nlon))
for it in range(0,nlat):
    for jt in range(0,nlon):
        if mask[it,jt] == 1.:
           coefa=np.zeros((3))
           coefb=np.zeros((3))
           hvar=np.zeros((3))
           tseries=cycle[:,it,jt]
           coefa,coefb,hvar=Harmonics(coefa,coefb,hvar,tseries,3,missval)
           harm1[it,jt]=hvar[0]
           harm2[it,jt]=hvar[1]
           harm3[it,jt]=hvar[2]
           harmonic1[:,it,jt]=rm[it,jt]
           harmonic1[:,it,jt]=harmonic1[:,it,jt]+coefa[0]*np.cos(2.*math.pi*time[:]/float(tot))+coefb[0]*np.sin(2.*math.pi*time[:]/float(tot))
           smoothed[:,it,jt]=np.mean(cycle[:,it,jt])
           for pp in range(0,3):
               smoothed[:,it,jt]=smoothed[:,it,jt]+coefa[pp]*np.cos(2.*math.pi*time[:]*(pp+1)/float(tot))+coefb[pp]*np.sin(2.*math.pi*time[:]*(pp+1)/float(tot))

print("Mean annual cycle smoothed.")
"""
Removing regions with zero or more than one rainy season per year

"""
id=np.where(harm2 >= harm1)
mask[id]=0.
rm[id]=0.
id=np.where(harm3 >= harm1)
mask[id]=0.
rm[id]=0.

"""
Calculating the day [day of year] that will be used as starting point (t0) for
the calculation of the rainy and dry seasons characteristics

"""
startwet=np.zeros((nlat,nlon))
for it in range(0,nlat):
    for jt in range(0,nlon):
        if mask[it,jt] == 1.:
           id=np.where(harmonic1[:,it,jt] == harmonic1[:,it,jt].min())
           startwet[it,jt]=jday[id]

#=======================================================================================
"""

Calculating the onset date of the rainy and dry seasons

"""
#=======================================================================================
cycle=None
smoothed=None
harmonic1=None

nyrs=int(year.max()-year.min())+1
ap=np.zeros((ntot))
onset_jday=np.zeros((nyrs,nlat,nlon))
demise_jday=np.zeros((nyrs,nlat,nlon))
onset_day=np.zeros((nyrs,nlat,nlon))
demise_day=np.zeros((nyrs,nlat,nlon))
onset_month=np.zeros((nyrs,nlat,nlon))
demise_month=np.zeros((nyrs,nlat,nlon))
onset_year=np.zeros((nyrs,nlat,nlon))
demise_year=np.zeros((nyrs,nlat,nlon))
durwet=np.zeros((nyrs,nlat,nlon))
durdry=np.zeros((nyrs,nlat,nlon))
totwet=np.zeros((nyrs,nlat,nlon))
totdry=np.zeros((nyrs,nlat,nlon))
wscurve=np.zeros((nyrs,int(tot/2),nlat,nlon))
dscurve=np.zeros((nyrs,int(tot/2),nlat,nlon))

wjd=np.zeros((nyrs))
wd=np.zeros((nyrs))
wm=np.zeros((nyrs))
wy=np.zeros((nyrs))
wsc=np.zeros((nyrs,int(tot/2)))

djd=np.zeros((nyrs))
dd=np.zeros((nyrs))
dm=np.zeros((nyrs))
dy=np.zeros((nyrs))
dsc=np.zeros((nyrs,int(tot/2)))

print('Calculating the onset of the rainy season...')

prec[prec<0.]=0.   # VERY IMPORTANT

npass=50
for it in range(0,nlat):
    print(it+1,' of ',nlat)
    for jt in range(0,nlon):
        if rm[it,jt] > 0.:
           sdate=startwet[it,jt]
           wjd[:]=0.
           wd[:]=0.
           wm[:]=0.
           wy[:]=0.
           wsc[:]=0.
           ap[:]=prec[:,it,jt]-rm[it,jt]
           wjd[:],wd[:],wm[:],wy[:],wsc[:,:]=rainyseason_onset(nyrs,tot,jday[:],day[:],month[:],year[:],sdate,ap[:],wjd[:],wd[:],wm[:],wy[:],wsc[:,:])
#------------------------------------------------------------------------
#    Quality control: removing outliers
#------------------------------------------------------------------------
# ---- Averagind dates using circular statistics
           miss=np.where(wjd[:] == 0.)
           id=np.where(wjd[:] != 0.)
           if len(id[0]) > 0:
              tmpx=np.cos(wjd[:]*math.pi/183.)
              tmpy=np.sin(wjd[:]*math.pi/183.)
              med=math.atan2(np.median(tmpy[id]),np.median(tmpx[id]))*183./math.pi
              if med < 0.:
                 med=med+tot
# ---- Converting dates close to beginning and end of the year
              tmpc=wjd[:]-med
              if len(miss[0]) > 0.:
                 tmpc[miss]=0.
              pos=np.where(tmpc[:] > float(tot)*0.5)
              if len(pos[0]) > 0:
                 tmpc[pos]=tmpc[pos]-tot
              neg=np.where(tmpc[:] < float(tot)*(-0.5))
              if len(neg[0]) > 0:
                 tmpc[neg]=tmpc[neg]+tot
# ---- Removing outliers (greater than 3 x IQR)
              iqr=np.percentile(tmpc[id],75)-np.percentile(tmpc[id],25)
              outl=np.where(np.abs(tmpc[:]) > iqr*1.5)
              if len(outl[0]) > 0:
                 wjd[outl]=0.
                 wd[outl]=0.
                 wm[outl]=0.
                 wy[outl]=0.
              onset_jday[:,it,jt]=wjd[:]
              onset_day[:,it,jt]=wd[:]
              onset_month[:,it,jt]=wm[:]
              onset_year[:,it,jt]=wy[:]
              wscurve[:,:,it,jt]=wsc[:,:]
#------------------------------------------------------------------------
#    Second pass (Bombardi et al. 2017)
#------------------------------------------------------------------------
              outl=np.where(wjd == 0.)
              if len(outl[0]) > 0:
                 wjd[:]=0.
                 wd[:]=0.
                 wm[:]=0.
                 wy[:]=0.
                 wsc[:]=0.
                 wjd[:],wd[:],wm[:],wy[:]=rainyseason_B17_onset(nyrs,tot,jday[:],day[:],month[:],year[:],sdate,ap[:],npass,wjd[:],wd[:],wm[:],wy[:])
#                 wjd[:],wd[:],wm[:],wy[:]=rainyseason_harmonic_onset(nyrs,tot,jday[:],day[:],month[:],year[:],sdate,ap[:],wjd[:],wd[:],wm[:],wy[:])
#------------------------------------------------------------------------
#    Quality control: removing outliers for the second pass
#------------------------------------------------------------------------
                 onset_jday[outl,it,jt]=wjd[outl]
                 onset_day[outl,it,jt]=wd[outl]
                 onset_month[outl,it,jt]=wm[outl]
                 onset_year[outl,it,jt]=wy[outl]
# ---- Averagind dates using circular statistics
                 wjd[:]=onset_jday[:,it,jt] #making sure to keep all the correct dates
                 miss=np.where(wjd[:] == 0.)
                 id=np.where(wjd[:] != 0.)
                 if len(id[0]) > 0:
                    tmpx=np.cos(wjd[:]*math.pi/183.)
                    tmpy=np.sin(wjd[:]*math.pi/183.)
                    med=math.atan2(np.median(tmpy[id]),np.median(tmpx[id]))*183./math.pi
                    if med < 0.:
                       med=med+tot
# ---- Converting dates close to beginning and end of the year
                    tmpc=wjd[:]-med
                    if len(miss[0]) > 0.:
                       tmpc[miss]=0.
                    pos=np.where(tmpc[:] > float(tot)*0.5)
                    if len(pos[0]) > 0:
                       tmpc[pos]=tmpc[pos]-tot
                    neg=np.where(tmpc[:] < float(tot)*(-0.5))
                    if len(neg[0]) > 0:
                       tmpc[neg]=tmpc[neg]+tot
# ---- Removing outliers (greater than 3 x IQR)
                    iqr=np.percentile(tmpc[id],75)-np.percentile(tmpc[id],25)
                    outl=np.where(np.abs(tmpc[:]) > iqr*3.)
                    if len(outl[0]) > 0:
                       onset_jday[outl,it,jt]=0.
                       onset_day[outl,it,jt]=0.
                       onset_month[outl,it,jt]=0.
                       onset_year[outl,it,jt]=0.
# print('Calculating the onset of the dry season...')
##----- Calculating the stats of the onset date of the dry season
#for it in range(0,nlat):
#    print(it+1,' of ',nlat)
#    for jt in range(0,nlon):
        if rm[it,jt] > 0.:
           sdate=startwet[it,jt] #It has to be wet because we are calculating it retrospectively
           djd[:]=0
           dd[:]=0
           dm[:]=0
           dy[:]=0
           dsc[:]=0
           ap[:]=prec[:,it,jt]-rm[it,jt]
#------------------------------------------------------------------------
#    First pass (Liebman & MArengo, 2001)
#------------------------------------------------------------------------
           djd[:],dd[:],dm[:],dy[:],dsc[:,:]=rainyseason_demise(nyrs,tot,jday[:],day[:],month[:],year[:],sdate,ap[:],djd[:],dd[:],dm[:],dy[:],dsc[:,:])
#------------------------------------------------------------------------
#    Quality control: removing outliers
#------------------------------------------------------------------------
# ---- Averagind dates using circular statistics
           miss=np.where(djd[:] == 0.)
           id=np.where(djd[:] != 0.)
           if len(id[0]) > 0:
              tmpx=np.cos(djd[:]*math.pi/183.)
              tmpy=np.sin(djd[:]*math.pi/183.)
              med=math.atan2(np.median(tmpy[id]),np.median(tmpx[id]))*183./math.pi
              if med < 0.:
                 med=med+tot
# ---- Converting dates close to beginning and end of the year
              tmpc=djd[:]-med
              if len(miss[0]) > 0.:
                 tmpc[miss]=0.
              pos=np.where(tmpc[:] > float(tot)*0.5)
              if len(pos[0]) > 0:
                 tmpc[pos]=tmpc[pos]-tot
              neg=np.where(tmpc[:] < float(tot)*(-0.5))
              if len(neg[0]) > 0:
                 tmpc[neg]=tmpc[neg]+tot
# ---- Removing outliers (greater than 3 x IQR)
              iqr=np.percentile(tmpc[id],75)-np.percentile(tmpc[id],25)
              outl=np.where(np.abs(tmpc[:]) > iqr*1.5)
              if len(outl[0]) > 0:
                 djd[outl]=0.
                 dd[outl]=0.
                 dm[outl]=0.
                 dy[outl]=0.
              demise_jday[:,it,jt]=djd[:]
              demise_day[:,it,jt]=dd[:]
              demise_month[:,it,jt]=dm[:]
              demise_year[:,it,jt]=dy[:]
              dscurve[:,:,it,jt]=dsc[:,:]
#------------------------------------------------------------------------
#    Second pass (Bombardi et al. 2017)
#------------------------------------------------------------------------
              outl=np.where(djd == 0.)
              if len(outl[0]) > 0:
                 djd[:]=0
                 dd[:]=0
                 dm[:]=0
                 dy[:]=0
                 dsc[:]=0
                 djd[:],dd[:],dm[:],dy[:]=rainyseason_B17_demise(nyrs,tot,jday[:],day[:],month[:],year[:],sdate,ap[:],npass,djd[:],dd[:],dm[:],dy[:])
#                 djd[:],dd[:],dm[:],dy[:]=rainyseason_harmonic_demise(nyrs,tot,jday[:],day[:],month[:],year[:],sdate,ap[:],djd[:],dd[:],dm[:],dy[:])
#------------------------------------------------------------------------
#    Quality control: removing outliers for the second pass
#------------------------------------------------------------------------
                 demise_jday[outl,it,jt]=djd[outl]
                 demise_day[outl,it,jt]=dd[outl]
                 demise_month[outl,it,jt]=dm[outl]
                 demise_year[outl,it,jt]=dy[outl]
# ---- Averagind dates using circular statistics
                 djd[:]=demise_jday[:,it,jt]
                 miss=np.where(djd[:] == 0.)
                 id=np.where(djd[:] != 0.)
                 if len(id[0]) > 0:
                    tmpx=np.cos(djd[:]*math.pi/183.)
                    tmpy=np.sin(djd[:]*math.pi/183.)
                    med=math.atan2(np.median(tmpy[id]),np.median(tmpx[id]))*183./math.pi
                    if med < 0.:
                       med=med+tot
# ---- Converting dates close to beginning and end of the year
                    tmpc=djd[:]-med
                    if len(miss[0]) > 0.:
                       tmpc[miss]=0.
                    pos=np.where(tmpc[:] > float(tot)*0.5)
                    if len(pos[0]) > 0:
                       tmpc[pos]=tmpc[pos]-tot
                    neg=np.where(tmpc[:] < float(tot)*(-0.5))
                    if len(neg[0]) > 0:
                       tmpc[neg]=tmpc[neg]+tot
# ---- Removing outliers (greater than 3 x IQR)
                    iqr=np.percentile(tmpc[id],75)-np.percentile(tmpc[id],25)
                    outl=np.where(np.abs(tmpc[:]) > iqr*3.)
                    if len(outl[0]) > 0:
                       demise_jday[outl,it,jt]=0.
                       demise_day[outl,it,jt]=0.
                       demise_month[outl,it,jt]=0.
                       demise_year[outl,it,jt]=0.
#------------------------------------------------------------------------
#    Masking regions where > 33% of the data are missing values 
#------------------------------------------------------------------------
           id=np.where(onset_jday[:,it,jt] == 0.)
           if len(id[0])/float(nyrs) > 0.33:
              rm[it,jt]=0.
              onset_jday[:,it,jt]=0.
              onset_day[:,it,jt]=0.
              onset_month[:,it,jt]=0.
              onset_year[:,it,jt]=0.
           id=np.where(demise_jday[:,it,jt] == 0.)
           if len(id[0])/float(nyrs) > 0.33:
              rm[it,jt]=0.
              demise_jday[:,it,jt]=0.
              demise_day[:,it,jt]=0.
              demise_month[:,it,jt]=0.
              demise_year[:,it,jt]=0.
# Rearranging years to account for retrospective calculation of demises
           if demise_year[1,it,jt] == float(yr0) or demise_year[2,it,jt] == float(yr0+1):
              demise_year[0:nyrs-1,it,jt]=demise_year[1:nyrs,it,jt]
              demise_year[nyrs-1,it,jt]=0.
              demise_month[0:nyrs-1,it,jt]=demise_month[1:nyrs,it,jt]
              demise_month[nyrs-1,it,jt]=0.
              demise_day[0:nyrs-1,it,jt]=demise_day[1:nyrs,it,jt]
              demise_day[nyrs-1,it,jt]=0.
              demise_jday[0:nyrs-1,it,jt]=demise_jday[1:nyrs,it,jt]
              demise_jday[nyrs-1,it,jt]=0.
#=======================================================================================
"""
Calculating duration of the wet and dry seasons and the total precipitated during the 
wet and dry seasons

"""
#=======================================================================================
           for yt in range(0,nyrs):
               if demise_year[yt,it,jt] == onset_year[yt,it,jt] and demise_year[yt,it,jt] != 0.:
                  if demise_jday[yt,it,jt] < onset_jday[yt,it,jt]:
                     #print(yt,"1 of 1")
#This means the dry season happens during the same year
                     beg=int((demise_year[yt,it,jt]-yr0)*tot+demise_jday[yt,it,jt]-1)
                     ned=int((onset_year[yt,it,jt]-yr0)*tot+onset_jday[yt,it,jt]-1)
                     durdry[yt,it,jt]=float(ned-beg)
# no (-1) because python doesn't use the last element anyway
                     totdry[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
                     # still have to calculate wet season properties
                     if yt < nyrs-1:
                        if demise_year[yt+1,it,jt] > 0.:
                           #print(yt,"3 of 1")
                           beg=int((onset_year[yt,it,jt]-yr0)*tot+onset_jday[yt,it,jt]-1)
                           ned=int((demise_year[yt+1,it,jt]-yr0)*tot+demise_jday[yt+1,it,jt]-1)
                           durwet[yt,it,jt]=float(ned-beg)
                           totwet[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
                  if onset_jday[yt,it,jt] < demise_jday[yt,it,jt]:
#This means the wet season happens during the same year
                     #print(yt,"1 of 1")
                     beg=int((onset_year[yt,it,jt]-yr0)*tot+onset_jday[yt,it,jt]-1)
                     ned=int((demise_year[yt,it,jt]-yr0)*tot+demise_jday[yt,it,jt]-1)
                     durwet[yt,it,jt]=float(ned-beg)
# no (-1) because python doesn't use the last element anyway
                     totwet[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
                     # still have to calculate dry season properties
                     if yt < nyrs-1:
                        if onset_year[yt+1,it,jt] > 0.:
                           #print(yt,"3 of 1")
                           beg=int((demise_year[yt,it,jt]-yr0)*tot+demise_jday[yt,it,jt]-1)
                           ned=int((onset_year[yt+1,it,jt]-yr0)*tot+onset_jday[yt+1,it,jt]-1)
                           durdry[yt,it,jt]=float(ned-beg)
                           totdry[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
               if 0. < demise_year[yt,it,jt] < onset_year[yt,it,jt]:
#this means the onset of the rainy season was found in year+1
                  #print(yt,"1 of 2")
                  beg=int((demise_year[yt,it,jt]-yr0)*tot+demise_jday[yt,it,jt]-1)
                  ned=int((onset_year[yt,it,jt]-yr0)*tot+onset_jday[yt,it,jt]-1)
                  durdry[yt,it,jt]=float(ned-beg)
                  totdry[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
                  if yt < nyrs-1:
                     if demise_year[yt+1,it,jt] > 0.:
                        if onset_jday[yt,it,jt] < demise_jday[yt+1,it,jt]:
                           #print(yt,"3 of 2")
                           beg=int((onset_year[yt,it,jt]-yr0)*tot+onset_jday[yt,it,jt]-1)
                           ned=int((demise_year[yt+1,it,jt]-yr0)*tot+demise_jday[yt+1,it,jt]-1)
                           durwet[yt,it,jt]=float(ned-beg)
                           totwet[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
               if 0. < onset_year[yt,it,jt] < demise_year[yt,it,jt]:
#this means the end of the rainy season was found in year+1
                  #print(yt,"1 of 3")
                  beg=int((onset_year[yt,it,jt]-yr0)*tot+onset_jday[yt,it,jt]-1)
                  ned=int((demise_year[yt,it,jt]-yr0)*tot+demise_jday[yt,it,jt]-1)
                  durwet[yt,it,jt]=float(ned-beg)
                  totwet[yt,it,jt]=np.sum(prec[beg:ned,it,jt])
                  if yt < nyrs-1:
                     if onset_year[yt+1,it,jt] > 0.:
                        if demise_jday[yt,it,jt] < onset_jday[yt+1,it,jt]:
                           #print(yt,"3 of 3")
                           beg=int((demise_year[yt,it,jt]-yr0)*tot+demise_jday[yt,it,jt]-1)
                           ned=int((onset_year[yt+1,it,jt]-yr0)*tot+onset_jday[yt+1,it,jt]-1)
                           durdry[yt,it,jt]=float(ned-beg)
                           totdry[yt,it,jt]=np.sum(prec[beg:ned,it,jt])

#=======================================================================================
"""
Preparing to save the results
"""
#=======================================================================================

print('Saving results...')
prec=None
dyr=np.arange(0,nyrs,1.)
dtime=np.arange(0,int(tot/2),1.)
yrs=np.arange(yr0,year[ntot-1]+1,1.)

onset_jday[onset_jday==0.]=missval
onset_day[onset_day==0.]=missval
onset_month[onset_month==0.]=missval
onset_year[onset_year==0.]=missval

demise_jday[demise_jday==0.]=missval
demise_day[demise_day==0.]=missval
demise_month[demise_month==0.]=missval
demise_year[demise_year==0.]=missval

totwet[totwet==0.]=missval
totdry[totdry==0.]=missval
durdry[durdry==0.]=missval
durwet[durwet==0.]=missval

#=======================================================================================
"""
Saving rainy and wet season characteristics
"""
#=======================================================================================

outfile = pathout+"onset.wet.season.CPC_UNI."+str(yrs[0])[0:4]+"-"+str(yrs[nyrs-3])[0:4]+".nc"
from netCDF4 import Dataset
rootgrp = Dataset(outfile, "w", format="NETCDF4")
rootgrp.close()
rootgrp = Dataset(outfile, "a")
# Creating dimensions
lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
lev = rootgrp.createDimension("lev", 1)
time = rootgrp.createDimension("time",nyrs-2)
#Creating coordinates
times = rootgrp.createVariable("time","i4",("time",))
levels = rootgrp.createVariable("lev","i4",("lev",))
latitudes = rootgrp.createVariable("lat","f8",("lat",))
longitudes = rootgrp.createVariable("lon","f8",("lon",))
# Filling coordinates
longitudes.units='degrees_east'
longitudes.long_name='Longitude'
longitudes[:]=lons[0:nlon]
latitudes.units='degrees_north'
latitudes.long_name='Latitude'
latitudes[:]=lats[0:nlat]
levels.units='millibar'
levels.long_name='Level'
levels[:]=1000.
times.long_name='Time'
times.units='years since '+str(yrs[0])[0:4]+'-01-01 00:00'
times[:]=dyr[0:nyrs-2]
# Creating variables 
var = rootgrp.createVariable("DOY","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var.long_name = 'Wet season onset [Day of Year]'
var[:,0,:,:] = onset_jday[0:nyrs-2,:,:]
# Creating variables 
var2 = rootgrp.createVariable("day","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var2.long_name = 'Wet season onset [Day of the month]'
var2[:,0,:,:] = onset_day[0:nyrs-2,:,:]
# Creating variables 
var3 = rootgrp.createVariable("month","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var3.long_name = 'Wet season onset month'
var3[:,0,:,:] = onset_month[0:nyrs-2,:,:]
# Creating variables 
var4 = rootgrp.createVariable("year","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var4.long_name = 'Wet season onset year'
var4[:,0,:,:] = onset_year[0:nyrs-2,:,:]


#=============
#  End of wet season/ Start of dry season
#=============

outfile = pathout+"demise.wet.season.CPC_UNI."+str(yrs[0])[0:4]+"-"+str(yrs[nyrs-3])[0:4]+".nc"
from netCDF4 import Dataset
rootgrp = Dataset(outfile, "w", format="NETCDF4")
rootgrp.close()
rootgrp = Dataset(outfile, "a")
# Creating dimensions
lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
lev = rootgrp.createDimension("lev", 1)
time = rootgrp.createDimension("time",nyrs-2)
#Creating coordinates
times = rootgrp.createVariable("time","i4",("time",))
levels = rootgrp.createVariable("lev","i4",("lev",))
latitudes = rootgrp.createVariable("lat","f8",("lat",))
longitudes = rootgrp.createVariable("lon","f8",("lon",))
# Filling coordinates
longitudes.units='degrees_east'
longitudes.long_name='Longitude'
longitudes[:]=lons[0:nlon]
latitudes.units='degrees_north'
latitudes.long_name='Latitude'
latitudes[:]=lats[0:nlat]
levels.units='millibar'
levels.long_name='Level'
levels[:]=1000.
times.long_name='Time'
times.units='years since '+str(yrs[0])[0:4]+'-01-01 00:00'
times[:]=dyr[0:nyrs-2]
# Creating variables 
var = rootgrp.createVariable("DOY","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var.long_name = 'Wet season demise [Day of Year]'
var[:,0,:,:] = demise_jday[0:nyrs-2,:,:]
# Creating variables 
var2 = rootgrp.createVariable("day","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var2.long_name = 'Wet season demise [Day of the month]'
var2[:,0,:,:] = demise_day[0:nyrs-2,:,:]
# Creating variables 
var3 = rootgrp.createVariable("month","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var3.long_name = 'Wet season demise month'
var3[:,0,:,:] = demise_month[0:nyrs-2,:,:]
# Creating variables 
var4 = rootgrp.createVariable("year","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var4.long_name = 'Wet season demise year'
var4[:,0,:,:] = demise_year[0:nyrs-2,:,:]

#=============
#  Total precipitation and duration
#=============

outfile = pathout+"total.precip.wet.season.CPC_UNI."+str(yrs[0])[0:4]+"-"+str(yrs[nyrs-2])[0:4]+".nc"
from netCDF4 import Dataset
rootgrp = Dataset(outfile, "w", format="NETCDF4")
rootgrp.close()
rootgrp = Dataset(outfile, "a")
# Creating dimensions
lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
lev = rootgrp.createDimension("lev", 1)
time = rootgrp.createDimension("time",nyrs)
#Creating coordinates
times = rootgrp.createVariable("time","i4",("time",))
levels = rootgrp.createVariable("lev","i4",("lev",))
latitudes = rootgrp.createVariable("lat","f8",("lat",))
longitudes = rootgrp.createVariable("lon","f8",("lon",))
# Filling coordinates
longitudes.units='degrees_east'
longitudes.long_name='Longitude'
longitudes[:]=lons[0:nlon]
latitudes.units='degrees_north'
latitudes.long_name='Latitude'
latitudes[:]=lats[0:nlat]
levels.units='millibar'
levels.long_name='Level'
levels[:]=1000.
times.long_name='Time'
times.units='years since '+str(yrs[0])[0:4]+'-01-01 00:00'
times[:]=dyr[0:nyrs]
# Creating variables 
var = rootgrp.createVariable("totwet","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var.long_name = 'Total precipiation during the wet season [mm]'
var[:,0,:,:] = totwet[0:nyrs,:,:]

outfile = pathout+"total.precip.dry.season.CPC_UNI."+str(yrs[0])[0:4]+"-"+str(yrs[nyrs-2])[0:4]+".nc"
from netCDF4 import Dataset
rootgrp = Dataset(outfile, "w", format="NETCDF4")
rootgrp.close()
rootgrp = Dataset(outfile, "a")
# Creating dimensions
lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
lev = rootgrp.createDimension("lev", 1)
time = rootgrp.createDimension("time",nyrs)
#Creating coordinates
times = rootgrp.createVariable("time","i4",("time",))
levels = rootgrp.createVariable("lev","i4",("lev",))
latitudes = rootgrp.createVariable("lat","f8",("lat",))
longitudes = rootgrp.createVariable("lon","f8",("lon",))
# Filling coordinates
longitudes.units='degrees_east'
longitudes.long_name='Longitude'
longitudes[:]=lons[0:nlon]
latitudes.units='degrees_north'
latitudes.long_name='Latitude'
latitudes[:]=lats[0:nlat]
levels.units='millibar'
levels.long_name='Level'
levels[:]=1000.
times.long_name='Time'
times.units='years since '+str(yrs[0])[0:4]+'-01-01 00:00'
times[:]=dyr[0:nyrs]
# Creating variables 
var = rootgrp.createVariable("totdry","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var.long_name = 'Total precipitation during the dry season [mm]'
var[:,0,:,:] = totdry[0:nyrs,:,:]


outfile = pathout+"duration.wet.season.CPC_UNI."+str(yrs[0])[0:4]+"-"+str(yrs[nyrs-2])[0:4]+".nc"
from netCDF4 import Dataset
rootgrp = Dataset(outfile, "w", format="NETCDF4")
rootgrp.close()
rootgrp = Dataset(outfile, "a")
# Creating dimensions
lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
lev = rootgrp.createDimension("lev", 1)
time = rootgrp.createDimension("time",nyrs)
#Creating coordinates
times = rootgrp.createVariable("time","i4",("time",))
levels = rootgrp.createVariable("lev","i4",("lev",))
latitudes = rootgrp.createVariable("lat","f8",("lat",))
longitudes = rootgrp.createVariable("lon","f8",("lon",))
# Filling coordinates
longitudes.units='degrees_east'
longitudes.long_name='Longitude'
longitudes[:]=lons[0:nlon]
latitudes.units='degrees_north'
latitudes.long_name='Latitude'
latitudes[:]=lats[0:nlat]
levels.units='millibar'
levels.long_name='Level'
levels[:]=1000.
times.long_name='Time'
times.units='years since '+str(yrs[0])[0:4]+'-01-01 00:00'
times[:]=dyr[0:nyrs]
# Creating variables 
var = rootgrp.createVariable("durwet","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var.long_name = 'Duration of the wet weson [day]'
var[:,0,:,:] = durwet[0:nyrs,:,:]


outfile = pathout+"duration.dry.season.CPC_UNI."+str(yrs[0])[0:4]+"-"+str(yrs[nyrs-2])[0:4]+".nc"
from netCDF4 import Dataset
rootgrp = Dataset(outfile, "w", format="NETCDF4")
rootgrp.close()
rootgrp = Dataset(outfile, "a")
# Creating dimensions
lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
lev = rootgrp.createDimension("lev", 1)
time = rootgrp.createDimension("time",nyrs)
#Creating coordinates
times = rootgrp.createVariable("time","i4",("time",))
levels = rootgrp.createVariable("lev","i4",("lev",))
latitudes = rootgrp.createVariable("lat","f8",("lat",))
longitudes = rootgrp.createVariable("lon","f8",("lon",))
# Filling coordinates
longitudes.units='degrees_east'
longitudes.long_name='Longitude'
longitudes[:]=lons[0:nlon]
latitudes.units='degrees_north'
latitudes.long_name='Latitude'
latitudes[:]=lats[0:nlat]
levels.units='millibar'
levels.long_name='Level'
levels[:]=1000.
times.long_name='Time'
times.units='years since '+str(yrs[0])[0:4]+'-01-01 00:00'
times[:]=dyr[0:nyrs]
# Creating variables 
var = rootgrp.createVariable("durdry","f4", ("time","lev","lat","lon",),fill_value=missval)
# Filling Variable
var.long_name = 'Duration of the dry season [day]'
var[:,0,:,:] = durdry[0:nyrs,:,:]

#========================================================================
#                             End of program
#========================================================================






