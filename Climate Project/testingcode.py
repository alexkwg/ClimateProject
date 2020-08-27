import netCDF4
import matplotlib.pyplot as plt
import datetime 
import statsmodels as statsmodels
import statistics as stats
import decimal
import os
import numpy as np
import copy

filename= 'TRMM1/3B42_Daily.19980101.7.nc4.nc4'
f=netCDF4.Dataset(filename)
lat, lon= [f.variables['lat'][:], f.variables['lon'][:]] #13 longitudes, 11 latitudes



class ReadData:
    """
    ReadData-- generates the dataset
    
    TODO: currently doesn't use generators, so it can be slow.
    """
    DMI_data = {}
    DMI_list = [] #sorted version of the dictionary
    sortedKeys_DMI = []
    TRMM_data = {}
    TRMM_list = []
    TRMM_mod7 = 0
    sortedKeys_TRMM = []
    ElNino_data= {}
    ElNino_list = [[],[],[],[]]
    sortedKeys_Nino = [[],[],[],[]] #currently, we have 4 Nino items
    data_dict = {}
    
    #constants
    begin = datetime.date(1998, 1, 1) #begin of dates we look at
    end = datetime.date(2019, 12, 30) #end of dates we look at
    
    
    
    
    def __init__(self):
        
        #in the init phase, we generated each data set and sort the dictionary keys
        
        self.DMI_data = self.generateDataDMI()
        print("DMI data generated")
        self.findMinData(self.DMI_data)
        self.sortedKeys_DMI = self.sortKeys(self.DMI_data)   
        for i in self.sortedKeys_DMI:
            self.DMI_list.append(self.DMI_data[i])
        
        self.TRMM_data = self.generateDataTrim()
        self.sortedKeys_TRMM = self.sortKeys(self.TRMM_data)  
        for i in self.sortedKeys_TRMM:
            self.TRMM_list.append(self.TRMM_data[i])
        print("TRMM data generated")
        sum_L = self.accumulateTRMM()
        mod7_0 = self.takeModN_TRMM(sum_L) #this creates a list with the summed precipitation fo the previous days
        #and we only take days that are mod7 == 0
        #
        
        
        
        self.ElNino_data = self.generateDataElNino()
        for i in range(0,len(self.ElNino_data)):
            self.sortedKeys_Nino[i] = self.sortKeys(self.ElNino_data[i])
            for j in self.sortedKeys_Nino[i]:
                self.ElNino_list[i].append(self.ElNino_data[i][j])
                #print(j)
        print("ELNino data generated")
        
    def accumulateTRMM(self):
        TRMM_list_update = []
        for i in range(6, len(self.TRMM_list)):
            TRMM_list_update.append(copy.deepcopy(self.TRMM_list[i]) ) #arrays are pass by REFRENCE--need to copy.
            for j in range(0,6):
                TRMM_list_update[i-6] += self.TRMM_list[i-j-1]
        return TRMM_list_update
    def takeModN_TRMM(self,L):
        TRMM_list_res = []
        for i in range(0,len(L)):
            if(i % 7 == self.TRMM_mod7): #only take number with mod of TRMM_mod7
                TRMM_list_res.append(L[i])
            
        return TRMM_list_res
    
        
    def sortKeys(self,data):
        return sorted(data,key = lambda key: self.toDateTime(key )-self.begin)
    def findMinData(self,data):
        v = 0
        for key in data.keys():
            diff = self.toDateTime(key )-self.begin
            if(v == 0 or v> diff):
                v = diff
        #print(v)
        return v
        
        
    
    """
    Generate DMI dataset--used for (???)
    """
    def generateDataDMI(self):
        iodfile = "dmi.nc" #were reading the DMI file
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

        date_DMI={}
        for j in range(0, len(SSTA)):
            if(self.isInTimeRange(DMIdate[j])):
                date_DMI[str(DMIdate[j]).replace("-","")] =  SSTA[j]
            else:
                pass
                #this data we will just discard for now
        return date_DMI
    
    def generateDataElNino(self):
        ninoindicesfile="elnino indices.txt"
        f=open(ninoindicesfile, "r")
        linelist=list(f.readlines())    
        #### ARRAYS OF FORMAT (DATETIME, DMI)
        date_Ninos=[{},{},{},{}]
        #12,3,34,4 are the ninos we take
        
        search_range = [[19,23],[32,36],[45,49],[58,62]]
        for j in range(0,len(search_range)):
            for i in range(4, len(linelist)):
                date_str = self.replaceDateTime(linelist[i][1:10])
                if (self.isInTimeRangeStr(date_str)):
                    date_Ninos[j][date_str] =  float(linelist[i][search_range[j][0]:search_range[j][1]])            
                else:
                    pass
                    #this is the date we discard
        f.close() 

        
        
        return date_Ninos
    """
    Generate TRMM dataset
    """
    def generateDataTrim(self): #this is very slow, potentially increase speed at some point
        TRMMfolder="TRMM1/"
        datetime_prec={} #create a dictionary with the dates
        for filename in os.listdir(TRMMfolder):
            if filename.endswith('nc4'):
                f=netCDF4.Dataset(TRMMfolder + str(filename))
                precp=f.variables['precipitation'][:]
                date_str = str(filename[11:19])
                #if(self.isInTimeRangeStr(date_str) and not str(self.toDateTime(date_str)-self.begin)=="0:00:00" and int(str(self.toDateTime(date_str)-self.begin).split(" ")[0]) % 7 == self.TRMM_mod7  ):
                if(self.isInTimeRangeStr(date_str)):    
                    datetime_prec[date_str] =  precp
                else:
                    pass
                f.close()
            else:
                pass
        
        #datetime_prec=np.array(datetime_prec)
        return datetime_prec
    def isInTimeRange(self,time):
        if (time >= self.begin and time <= self.end):
            return True
        return False
    def toDateTime(self,time):
        return datetime.date(int(time[0:4]),int(time[4:6]),int(time[6:8]))
    def isInTimeRangeStr(self,time):
        return self.isInTimeRange(self.toDateTime(time))
        
        
    def replaceDateTime(self,string):
        return str(datetime.datetime.strptime(string, '%d%b%Y')).split(" ")[0].replace("-","")
    




reader = ReadData()


class AnalyzeData:
    "class for anaylzing the data we generated from a reader"
    reader_class = None
    snhtsetDMI = []
    snhtsetNino = []
    def __init__(self,reader):
        self.reader_class = reader
        self.generateAllSNHT()
    def generateAllSNHT(self):
        self.snhtsetDMI = self.SNHT(reader.DMI_list)
        #snhtsetTRMM = self.SNHT(reader.TRMM_list) 
        #^^ we can't use this set, since we have more than one varaible. 
        #QUESTION: how should we deal with this set?   ## Run through a "multidimensional test". 
        self.snhtsetNino = []
        for i in reader.ElNino_list:    
            self.snhtsetNino.append(self.SNHT(i))
        #TODO: continue with taking the snhtset
    def SNHT(self,dataset):   #SNHT FOR time series of dimension 1 x n
        #This function is fairly slow, but we only need to run it once, so speed isn't all that important
        n=len(dataset)  
        sd=decimal.Decimal(stats.pstdev(dataset))
        mean = decimal.Decimal(stats.mean(dataset))
        snhtset=[]
        for y in range(1,n-1): #Loops from Week 1 to week n-1
            summ1=0 
            for i in range(1, y+1): #QUESTION: do we want to skip the first week? Why not include 0?
                summ1 += decimal.Decimal((decimal.Decimal(dataset[i]) - mean)/sd)
            z_1= decimal.Decimal(1/(y+1)) * summ1
            summ2=0
            for i in range(y+1, n):
                summ2 += decimal.Decimal((decimal.Decimal(dataset[i]) - mean)/sd)
            z_2=decimal.Decimal(1.0/(n-y)) *summ2
            snhtset.append(y * (z_1 **2) + (n-y)*(z_2**2))
        return snhtset

def Pearlson(x, y, lag):      #LAG means that we compare x_t against y_{t-lag} where lag is positive
    x=x[lag:]
    y=y[:-lag]
    x_mean=stats.mean(x)
    y_mean=stats.mean(y)
    top=0
    n=len(x)-lag
    for j in range(0, n):
        top=top+(x[j]-x_mean)*(y[j]-y_mean)
    pearlson= top/(n*stats.stdev(x) * stats.stdev(y))   
    return pearlson
a=reader.DMI_list
b=reader.ElNino_list[0]
t=reader.DMI_data.keys()
plt.plot(t, a, 'r', t, b, 'g')
plt.show()

