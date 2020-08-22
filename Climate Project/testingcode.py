import netCDF4
import matplotlib.pyplot as plt
import datetime 
import statsmodels as statsmodels
import statistics as stats
import xarray as xr
import pandas as pd
import decimal
import os
import numpy as np



######## DMI ARRAY CREATION #######
######## DMI ARRAY CREATION #######
######## DMI ARRAY CREATION #######
######## DMI ARRAY CREATION #######

iodfile = "dmi.nc"
f = netCDF4.Dataset(iodfile, 'r')
SSTAnp=f.variables['DMI'][:]
weeknp=f.variables['WEDCEN2'][:]
SSTA=SSTAnp.tolist() #list of DMI SSTA
week=weeknp.tolist() 
f.close()
DMIdate=[]     #DATE OF DMI 
for i in range(0, len(week)):    
    days=week[i]
    start=datetime.date(1900, 1, 1)
    delta=datetime.timedelta(days)
    offset=start+delta
    DMIdate.append(offset)


date_DMI=[]
for j in range(0, len(SSTA)):
    date_DMI.append([DMIdate[j], SSTA[j]])
date_DMI=np.array([date_DMI])
#print(date_DMI)    #### ARRAY OF FORMAT (DATETIME, DMI)


###### TRMM DATA MEGA ARRAY CREATION ######
###### TRMM DATA MEGA ARRAY CREATION ######
###### TRMM DATA MEGA ARRAY CREATION ######
###### TRMM DATA MEGA ARRAY CREATION ######
###### TRMM DATA MEGA ARRAY CREATION ######

TRMMfolder="TRMM1/"
datetime_prec=[]
for filename in os.listdir(TRMMfolder):
    if filename.endswith('nc4'):
        f=netCDF4.Dataset(TRMMfolder + str(filename))
        precp=f.variables['precipitation'][:]
        datetime_prec.append([filename[11:19], precp])
        f.close()
    else:
        pass
datetime_prec=np.array(datetime_prec)
print("Done with array creation")
print(datetime_prec[0][1][0] )





##### NINO ARRAY CREATION ######
##### NINO ARRAY CREATION ######
##### NINO ARRAY CREATION ######
##### NINO ARRAY CREATION ######
ninoindicesfile="elnino indices.txt"

f=open(ninoindicesfile, "r")
linelist=list(f.readlines())    
date_Nino12=[]                     #### ARRAY OF FORMAT (DATETIME, DMI)
for i in range(4, len(linelist)):
    date_Nino12.append([datetime.datetime.strptime(linelist[i][1:10], '%d%b%Y'), decimal.Decimal(linelist[i][19:23])])
    
    
date_Nino3=[]                     #### ARRAY OF FORMAT (DATETIME, DMI)
for i in range(4, len(linelist)):
    date_Nino3.append([datetime.datetime.strptime(linelist[i][1:10], '%d%b%Y'), decimal.Decimal(linelist[i][32:36])])
    
    
date_Nino34=[]                     #### ARRAY OF FORMAT (DATETIME, DMI)
for i in range(4, len(linelist)):
    date_Nino34.append([datetime.datetime.strptime(linelist[i][1:10], '%d%b%Y'), decimal.Decimal(linelist[i][45:49])])
    
    
date_Nino4=[]                     #### ARRAY OF FORMAT (DATETIME, DMI)
for i in range(4, len(linelist)):
    date_Nino4.append([datetime.datetime.strptime(linelist[i][1:10], '%d%b%Y'), decimal.Decimal(linelist[i][58:62])])

f.close()    
##### Defining STATISTICAL TESTS #####
##### Defining STATISTICAL TESTS #####
##### Defining STATISTICAL TESTS #####
##### Defining STATISTICAL TESTS #####    
##### Defining STATISTICAL TESTS #####
def SNHT(dataset):   #SNHT FOR time series of dimension 1 x n
    n=len(dataset)
    sd=stats.pstdev(dataset)
    snhtset=[]
    for y in range(1,n-1): #Loops from Week 1 to week n-1
        summ1=0
        for i in range(1, y+1):
            summ1 += decimal.Decimal((dataset[i] - decimal.Decimal(stats.mean(dataset)))/sd)
        z_1= decimal.Decimal(1/(y+1)) * summ1
        summ2=0
        for i in range(y+1, n):
            summ2 += decimal.Decimal((dataset[i] - decimal.Decimal(stats.mean(dataset)))/sd)
        z_2=decimal.Decimal(1/(n-y)) *summ2
        snhtset.append(y * (z_1 **2) + (n-y)*(z_2**2))
    return snhtset

######### GRAPHS###########
######### GRAPHS###########
######### GRAPHS###########

#plt.plot(DMIdate,SSTA)
#plt.show()

