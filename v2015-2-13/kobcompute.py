import re
import sys
import math
import koblib
import vecmath as vec

def upper(s):
    multipleofsmallest = s
    done = False
    counter = 0
    temp = multipleofsmallest
    while not done and s != 0:
        if temp % 10 != 0:
            done = True
        else:
            temp = temp / 10
            counter+=1 
    return multipleofsmallest - temp % 10 * (10 ** counter) + 10 ** (counter + 1)

def kobclean(statesL,currenttime):
    for state in statesL[1:]:
            s = state[0]
            if upper(s) == currenttime:
                statesL.remove(state)                                     

def countatoms(map, atomIDlist):
    sum=0
    for atomID in atomIDlist:
        sum+=len(map[atomID])    
    return sum       

def out_C_sq(filename, starttime):
    f = open(filename,'rU')
    timeskipto(f,starttime)
    fout = open('c_sq.txt','w')
 
    timestep = starttime

    for line in f:
        match0 = re.match("TIMESTEP (\d+)",line)
        match = re.match("EIGEN_VECTOR",line)
        if match0:
            timestep = match0.group(1)
 
        if match and int(timestep) % 100 == 0:
           L = []
           for c in f.next().split():
             L.append(float(c) * float(c))
           L = sorted(L, reverse=True)
           fout.write("%f   %f \n" % (L[0], L[1]))
           

 
def minimagedist(r2,r1,boxlength):
    diff = r2 - r1
    if math.fabs(diff) > boxlength / 2:
        if diff > 0:
            return diff - boxlength
        else:
            return diff + boxlength
    else:
        return diff

def minimage3Ddist(r2,r1,boxlength):
    sumsq=0
    for i in range(3):
        diff=r2[i]-r1[i]
        if math.fabs(diff) > boxlength/2:
            if diff>0:
                diff=diff-boxlength
            else:
                diff=diff+boxlength
        sumsq+=diff*diff
    return math.sqrt(sumsq) 

def calselfsqdist(statesL,previndex,currindex,typeid,q):
    Typemap = q.Typemap
    noofatom = statesL[currindex][1]    
    sqsum = 0
    
    timediff = statesL[currindex][0] - statesL[previndex][0]
    if q.correctdrift == False:
        for i in Typemap[typeid]:
              xdiff = statesL[currindex][i][2] - statesL[previndex][i][2] + (statesL[currindex][i][5] - statesL[previndex][i][5]) * q.xboxlength
              ydiff = statesL[currindex][i][3] - statesL[previndex][i][3] + (statesL[currindex][i][6] - statesL[previndex][i][6]) * q.yboxlength
              zdiff = statesL[currindex][i][4] - statesL[previndex][i][4] + (statesL[currindex][i][7] - statesL[previndex][i][7]) * q.zboxlength
              sqsum+=xdiff ** 2 + ydiff ** 2 + zdiff ** 2    
    else:
        for i in Typemap[typeid]:
              xdiff = statesL[currindex][i][2] - statesL[previndex][i][2] + (statesL[currindex][i][5] - statesL[previndex][i][5]) * q.xboxlength - q.driftvL[0] * timediff
              ydiff = statesL[currindex][i][3] - statesL[previndex][i][3] + (statesL[currindex][i][6] - statesL[previndex][i][6]) * q.yboxlength - q.driftvL[1] * timediff
              zdiff = statesL[currindex][i][4] - statesL[previndex][i][4] + (statesL[currindex][i][7] - statesL[previndex][i][7]) * q.zboxlength - q.driftvL[2] * timediff
              sqsum+=xdiff ** 2 + ydiff ** 2 + zdiff ** 2
    
    if sqsum > len(Typemap[typeid]) * 6 * 0.005 * timediff:
          print("current: %d previous: %d  , No. of Atoms of this type: %d , sqsum: %f" % (statesL[currindex][0],statesL[previndex][0],len(Typemap[typeid]),sqsum))                 
              
    return (sqsum,len(Typemap[typeid]))

def findwhichtimeslot(timediff,smallesttime):
    if timediff == 0:
        return -1
    else:
        time_unit = timediff / smallesttime #e.g timediff=400000, smallesttime =1000, time_unit=400
        whichpow = math.floor(math.log10(time_unit)) #whichpow=2
        firstdigit = time_unit / 10 ** whichpow #firstdigit=4
        return int(whichpow * 10 + firstdigit - (whichpow + 1)) #return 23

def calclastMSD(q,type):
    sum = 0
    sqdistinfo = calselfsqdist(q.statesL,0,-1,type,q)
    timediff = q.statesL[-1][0] - q.statesL[0][0]
    return (timediff,sqdistinfo[0] * 1. / sqdistinfo[1])


def calcMSD(q,MSDL,type):   
    no_of_timeslots = len(q.statesL)    
    currtimediff = q.statesL[-1][0] - q.statesL[0][0]    
    if no_of_timeslots > 0:
         for i in range(no_of_timeslots - 2,-1,-1): #reading the times backards
            timediff = q.statesL[-1][0] - q.statesL[i][0]
            print("%d:%d %d:%d   %d" % (no_of_timeslots - 1,q.statesL[-1][0],i, q.statesL[i][0],timediff))            
                    
            if timediff > 0:
                basetime = int(10 ** math.floor(math.log10(timediff)))
            else:
                sys.exit("Error: time difference >= 0") 
                
            
            if timediff % basetime == 0:           
              sqdistinfo = calselfsqdist(q.statesL,i,-1,type,q)
              if timediff == currtimediff:                  
                  MSDL.append([timediff,sqdistinfo[0],sqdistinfo[1]]) #sqdistinfo[0]=sqdist, [1]=count
              elif timediff < currtimediff:
                  slot = findwhichtimeslot(timediff,q.smalleststep)                
                  MSDL[slot + 1][1]+=sqdistinfo[0] #the +1 is for the label
                  MSDL[slot + 1][2]+=sqdistinfo[1]                        
            else:
                print timediff
               
     
def calc1DSqdisp(q,sqdispL,type):   
    no_of_timeslots = len(q.statesL)    
    currtimediff = q.statesL[-1][0] - q.statesL[0][0]    
    if no_of_timeslots > 0:
         for i in range(no_of_timeslots - 2,-1,-1): #reading the times backards
            timediff = q.statesL[-1][0] - q.statesL[i][0]
            print("%d:%d %d:%d   %d" % (no_of_timeslots - 1,q.statesL[-1][0],i, q.statesL[i][0],timediff))            
                    
            if timediff > 0:
                basetime = int(10 ** math.floor(math.log10(timediff)))
            else:
                sys.exit("Error: time difference >= 0") 

            if timediff % basetime == 0:                
              dispL = None
              slot = findwhichtimeslot(timediff,q.smalleststep) 
              if timediff == currtimediff:
                   dispL=[0]
                   sqdispL.append([timediff,dispL])
                                         
              elif timediff < currtimediff:                   
                   dispL=sqdispL[slot+1][1]
              
              currindex = -1
              previndex = i              
              for k in q.Typemap[type]:
                 xdiff = q.statesL[currindex][k][2] - q.statesL[previndex][k][2] + (q.statesL[currindex][k][5] - q.statesL[previndex][k][5]) * q.xboxlength                                  
                 ydiff = q.statesL[currindex][k][3] - q.statesL[previndex][k][3] + (q.statesL[currindex][k][6] - q.statesL[previndex][k][6]) * q.yboxlength
                 zdiff = q.statesL[currindex][k][4] - q.statesL[previndex][k][4] + (q.statesL[currindex][k][7] - q.statesL[previndex][k][7]) * q.zboxlength
                 whichbin=int(math.fabs(xdiff)/q.binsize)+1
                 koblib.histo(sqdispL[slot+1][1],q.binsize,whichbin,1)
                 dispL[0]+=1
                 whichbin=int(math.fabs(ydiff)/q.binsize)+1
                 koblib.histo(sqdispL[slot+1][1],q.binsize,whichbin,1)
                 dispL[0]+=1
                 whichbin=int(math.fabs(zdiff)/q.binsize)+1
                 koblib.histo(sqdispL[slot+1][1],q.binsize,whichbin,1)               
                 dispL[0]+=1
                 
def calselfangle(statesL,prev2index, previndex,currindex,typeid,q,dataL):
    Typemap_prev2 = statesL[prev2index][2][typeid]
    Typemap_prev = statesL[previndex][2][typeid]
    Typemap_curr = statesL[currindex][2][typeid]
    binsize=0.02
    
    noofele=len(Typemap_curr)    
    xl=q.xl
    
    for count in range(noofele):
        i=Typemap_curr[count]
        j=Typemap_prev[count]
        k=Typemap_prev2[count]
        pointi=[]
        for h in range(3):
            pointi.append(statesL[currindex][i][xl+h] +statesL[currindex][i][xl+3+h]*q.xboxlength)
        pointj=[]
        for h in range(3):
            pointj.append(statesL[previndex][j][xl+h] +statesL[previndex][j][xl+3+h]*q.xboxlength)
        pointk=[]
        for h in range(3):
            pointk.append(statesL[prev2index][k][xl+h] +statesL[prev2index][k][xl+3+h]*q.xboxlength)
        v1=[]
        v2=[]
        for h in range(3):
            v1.append(pointi[h]-pointj[h])
            v2.append(pointj[h]-pointk[h])        
        angle=vec.anglev1v2(v1,v2)
        whichbin=int(angle/3.141592654/2/binsize)+1 # +1 to account for the label slot
        koblib.histo(dataL,binsize,whichbin,1)  
    
def findlastangletime(statesL,startindex,time):      
    for i in range(startindex-1,-1,-1):
        if statesL[i][0]==time:
            return i
        elif statesL[i][0]<time:
            return -99
        

def calcangle(q,angleL,type):   
    no_of_timeslots = len(q.statesL)    
    currtimediff = q.statesL[-1][0] - q.statesL[0][0]    
    if no_of_timeslots > 0:
         for previndex in range(no_of_timeslots - 2,-1,-1): #reading the times backards
            timediff = q.statesL[-1][0] - q.statesL[previndex][0]
            #print("%d:%d %d:%d   %d" % (no_of_timeslots - 1,q.statesL[-1][0],i, q.statesL[i][0],timediff))            
             
            if timediff > 0 and timediff%1==0:
                basetime = int(10 ** math.floor(math.log10(timediff)))            
                if timediff % basetime == 0:                              
                   searchtime=q.statesL[previndex][0]-timediff
                   if searchtime>q.statesL[0][0]:
                       prev2index=findlastangletime(q.statesL,previndex,searchtime)
                       if prev2index!=-99 and timediff in q.angletimeL:       
                          listindex=q.angletimeL.index(timediff)
                          calselfangle(q.statesL,prev2index, previndex,-1,type,q,angleL[listindex])
                  
                               


                                         

def calcRDF(q,RDFL):
    pair = RDFL[0]
    list1=list2=None
    if len(pair[0])==1:
        list1 = q.Typemap[pair[0][0]]
    else:
        list1=[]
        for index in pair[0][:]:
            for item in q.Typemap[index][:]:
                list1.append(item)

    if len(pair[1])==1:
        list2 = q.Typemap[pair[1][0]]
    else:
        list2=[]
        for index in pair[1][:]:
            for item in q.Typemap[index][:]:
                list2.append(item)

    for i in list1:
            #print i
            atominfo_i=q.statesL[-1][i]
            for j in list2:              
              if j != i:                         
    
                xdiff = minimagedist(q.statesL[-1][j][2],atominfo_i[2],q.xboxlength)
                ydiff = minimagedist(q.statesL[-1][j][3],atominfo_i[3],q.yboxlength)
                zdiff = minimagedist(q.statesL[-1][j][4],atominfo_i[4],q.zboxlength)
                dist = math.sqrt(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)
                
                              
                whichbin = int(dist / q.binsize) + 1 # +1 to account for the label slot
                koblib.histo(RDFL,q.binsize,whichbin,1)
   

def calcClosest(q,ClosestL,listcount):
    pair = ClosestL[0]
    list1 = q.Typemap[pair[0][0]]    
    for item in pair[0][1:]:
        list1+=q.Typemap[item]

    list2 = []
    if pair[1][0] != 999:
        list2 = q.Typemap[pair[1][0]]
        for item in pair[1][1:]:
            list2+=q.Typemap[item]
    else:        
        for i in range(1,len(q.Typemap) - 1):
            if i not in notpolymerL:
                list2 = list2 + q.Typemap[i]         
        if needCEC == True: #if CEC, the last atom type is CEC and not polymer
            list2 = list2 + q.Typemap[-1]         
        elif (i + 1) not in notpolymerL:
                list2 = list2 + q.Typemap[-1]            
            
        
    
    closestindex = -1
    for i in list1:            
            mindist = q.zboxlength            
            for j in list2:
                 if j != i:                        
                       xdiff = minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],q.xboxlength)
                       ydiff = minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],q.yboxlength)
                       zdiff = minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],q.zboxlength)
                       dist = math.sqrt(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)
                       if dist < mindist:
                              mindist = dist                                                 
            whichbin = int(mindist / q.binsize) + 1 # +1 to account for the label slot
            koblib.histo(ClosestL,q.binsize,whichbin,1)
            if pair[2] == True:                
                list3 = q.Typemap[q.conditionL[listcount - 1][1][0]]
                cutoff = q.conditionL[listcount - 1][2]
                atomcount = 0
                for k in list3:
                    if k != i:
                       xdiff = minimagedist(q.statesL[-1][i][2],q.statesL[-1][k][2],q.xboxlength)
                       ydiff = minimagedist(q.statesL[-1][i][3],q.statesL[-1][k][3],q.yboxlength)
                       zdiff = minimagedist(q.statesL[-1][i][4],q.statesL[-1][k][4],q.zboxlength)
                       dist = math.sqrt(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)
                       if dist < cutoff:
                           atomcount+=1

                binL = q.bigConditionalL[listcount]
                maxbin = len(binL) - 1 
                if whichbin > maxbin:                      
                        while maxbin < whichbin: # -1 is for the label for the pair at the first slot
                            if len(binL) == 1:
                              binL.append([q.binsize / 2,0,0])
                            else:
                              binL.append([binL[-1][0] + q.binsize,0,0])  
                            maxbin+=1
                        binL[-1][1]+=atomcount #last entry
                        binL[-1][2]+=1 #last entry
                else:
                      binL[whichbin][1]+=atomcount    
                      binL[whichbin][2]+=1

        

count1 = 0
count2 = 0

def calcTseContourPlot(q,L1,L2,L3):
            
    count = 0
    for i in L1:
        
        j = L2[count]
        k = L3[count]
        
        #print q.statesL[-1][j][0],q.statesL[-1][k][0],q.statesL[-1][i][0]
        
        xdiff1 = minimagedist(q.statesL[-1][i][2],q.statesL[-1][j][2],q.xboxlength)
        ydiff1 = minimagedist(q.statesL[-1][i][3],q.statesL[-1][j][3],q.yboxlength)
        zdiff1 = minimagedist(q.statesL[-1][i][4],q.statesL[-1][j][4],q.zboxlength)
        dist = math.sqrt(xdiff1 ** 2 + ydiff1 ** 2 + zdiff1 ** 2)
        
        xdiff2 = minimagedist(q.statesL[-1][k][2],q.statesL[-1][j][2],q.xboxlength)
        ydiff2 = minimagedist(q.statesL[-1][k][3],q.statesL[-1][j][3],q.yboxlength)
        zdiff2 = minimagedist(q.statesL[-1][k][4],q.statesL[-1][j][4],q.zboxlength)
        
        temp1 = [xdiff1,ydiff1,zdiff1]
        len1 = koblib.lenvec(temp1)
        SC1vec = [xdiff1 / len1,ydiff1 / len1,zdiff1 / len1]
        temp2 = [xdiff2,ydiff2,zdiff2]
        len2 = koblib.lenvec(temp2)
        Cb1C2vec = [xdiff2 / len2,ydiff2 / len2,zdiff2 / len2]
       
        whichbin = int(dist / q.binsize) + 1 # +1 to account for the label slot
        
        angle = math.acos(koblib.dotvec(SC1vec,Cb1C2vec))
        whichanglebin = int(angle / q.anglebinsize)
        maxbin = len(q.C1SL) - 1
        if whichanglebin < 80:
          if whichbin > maxbin:
            while maxbin < whichbin: # -1 is for the label for the pair at the first slot
                if len(q.C1SL) == 1:
                    q.C1SL.append([q.binsize / 2,0])
                      
                else:
                    q.C1SL.append([q.C1SL[-1][0] + q.binsize,0])                     
                    maxbin+=1
                for i in range(0,80):
                    q.C1SL[-1].append(0)
            q.C1SL[-1][1]+=1 #last entry
            q.C1SL[-1][1 + 1 + whichanglebin]+=1 #last entry
          else:
            q.C1SL[whichbin][1]+=1
            q.C1SL[whichbin][1 + 1 + whichanglebin]+=1
        else:
           print whichanglebin
           print "Not collected"
        count+=1

def collect_ci(ciL,q,noofstates):
    
    q.bigciL[0]+=1
    for i in range(2):              
        cisq = ciL[i] ** 2
        whichbin = int(cisq / q.ci_binsize) + 1 # +1 to account for the label slot
        koblib.histo(q.bigciL[i + 1],q.ci_binsize,whichbin,1)         
    ci_sqdiff = ciL[0] * ciL[0] - ciL[1] * ciL[1]
    whichbin = int(ci_sqdiff/q.ci_binsize) +1 # +1 to account for the label slot
    koblib.histo(q.bigciL[3],q.ci_binsize,whichbin,1)  
    whichbin = int(ci_sqdiff / q.ci_binsize) + 1 # +1 to account for the label slot
    whichbin = noofstates + 1  # +1 to account for the label slot   
    koblib.histo(q.NstateL,1,whichbin,1)

def centralboxdist(q,atomstate1,atomstate2):
    xdiffsq=   atomstate1[2] + (atomstate1[q.xl+3] - atomstate2[q.xl+3]) * q.xboxlength
    ydiffsq=   atomstate1[3] + (atomstate1[q.xl+4] - atomstate2[q.xl+4]) * q.yboxlength
    zdiffsq=   atomstate1[4] + (atomstate1[q.xl+5] - atomstate2[q.xl+5]) * q.zboxlength   
   
    return math.sqrt(xdiffsq+ydiffsq+zdiffsq)

def endtoenddist(q,currentdist):
    #these indexes 3,125,247...are the head nad tail indexes of the chains.
    count=0
    sum=0
    for index in q.headtype:
        index2=q.tailtype[count]
        dist=centralboxdist(q,q.statesL[-1][index],q.statesL[-1][index2])        
        q.currentdist[count][0]+=dist
        q.currentdist[count][1]+=dist**2
        q.currentdist[count][2]+=1
        count+=1    
    return currentdist

def findwaterneighbors(seed,waterL,radius,q,localclusterL):
    seedwaterindex=waterL[seed]
    seedpos=q.statesL[-1][seedwaterindex][q.xl:q.xl+3] 
    ID=q.statesL[-1][seedwaterindex][0]
    localclusterL.append(q.statesL[-1][seedwaterindex][0])
    waterL[seed]=-999        
    for i in range(len(waterL)):
        #print i
        if waterL[i]!=-999:            
            index=waterL[i]  
            otherpos=q.statesL[-1][index][q.xl:q.xl+3]
            dist=minimage3Ddist(seedpos,otherpos,q.boxlength)
            if dist<radius:                    
                ID2=q.statesL[-1][index][0]                            
                findwaterneighbors(i,waterL,radius,q,localclusterL)

def findwaterclusters(q):
    radius=3.5
    waterL=q.Typemap[q.OWID][:]
    noofclusters=0     
    totalnumber=0   
    while len(waterL)>0:
        seedwaterindex=0        
        tempL=[]
        localclusterL=[]
        findwaterneighbors(seedwaterindex,waterL,radius,q,localclusterL)                           
        for item in waterL:
            if item!=-999:
                tempL.append(item)
        waterL=tempL
        clustersize=len(localclusterL)
        print clustersize
        totalnumber+=clustersize
        koblib.histo(q.watersizeclustersL,1,clustersize,1)
        noofclusters+=1
    koblib.histo(q.nowaterclustersL,1,noofclusters,1)
    print totalnumber

def findFPT(q):
    count=0
    currtime=q.statesL[-1][0]
    firstt=q.statesL[0][0]    
      
    for item in q.FPTtypes:
        MFPT=item[2]
        if MFPT!=-999:
            q.timeseries[count].append([currtime,[]])
        type=item[0]
        threshold=item[1]
        if len(q.bigtrackingFPTL[count])!=len(q.Typemap[type])+1: #+1 for label
            for item in q.Typemap[type]:
                q.bigtrackingFPTL[count].append([currtime,q.statesL[-1][item][q.xl:q.xl+3]])            
        else:
            index=1
            for item in q.Typemap[type]:
                currpos=q.statesL[-1][item][q.xl:q.xl+3]
                oldtime=q.bigtrackingFPTL[count][index][0]
                oldpos=q.bigtrackingFPTL[count][index][1]               
                dist=minimage3Ddist(currpos,oldpos,q.boxlength)
                if dist>=threshold:
                    timediff=currtime-oldtime
                    whichbin=timediff/q.FPTbinsize+1
                    koblib.histo(q.bigFPTL[count+1],q.FPTbinsize,whichbin,1)
                    q.bigtrackingFPTL[count][index][0]=currtime
                    q.bigtrackingFPTL[count][index][1]=currpos
                    if MFPT!=-999:
                        if timediff>=MFPT:
                            whichtimeslot=(oldtime-firstt)/q.smalleststep-1
                            howmanyslots=timediff/q.smalleststep+1
                            for i in range(howmanyslots):
                                q.timeseries[count][whichtimeslot+i][1].append(item)
                index+=1
        count+=1



def compute(q):         
    if q.MSDdecomposition:
        noofatom = q.statesL[-1][1]        
        distL = []
        contL = []
        if len(q.statesL[0]) == 3 + noofatom + q.noofCEC:
            for i in range(0,q.noofCEC):
                print i
                CECindex = noofatom + 3 + i
                q.statesL[0].append([q.statesL[-1][CECindex][0],-2,0,0,0,0,0,0])
            for i in range(0,q.noofCEC):
                CECindex = noofatom + 3 + i
                q.statesL[0].append([q.statesL[-1][CECindex][0],-3,0,0,0,0,0,0])
        for i in range(0,q.noofCEC):
          CECindex = noofatom + 3 + i
          xdiff = q.statesL[-1][CECindex][2] - q.statesL[-2][CECindex][2] + (q.statesL[-1][CECindex][5] - q.statesL[-2][CECindex][5]) * q.xboxlength
          ydiff = q.statesL[-1][CECindex][3] - q.statesL[-2][CECindex][3] + (q.statesL[-1][CECindex][6] - q.statesL[-2][CECindex][6]) * q.yboxlength
          zdiff = q.statesL[-1][CECindex][4] - q.statesL[-2][CECindex][4] + (q.statesL[-1][CECindex][7] - q.statesL[-2][CECindex][7]) * q.zboxlength

          reacted = False
          CEC2ndorderL = sorted(q.CECorderL, key=lambda k:k[1])
          #if q.CECorderL[i][2]==1:
          if CEC2ndorderL[i][2] == 1:
              reacted = True
          q.hopsL[i].append(q.statesL[-1][CECindex][0])
          if reacted == True:
              if q.statesL[-1][0] - q.statesL[0][0] < 500000:
                 q.counthopsL[i][0]+=1
                 q.counthopsL[i][2].append(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)
              distL.append([q.statesL[-1][CECindex][0],-3,q.statesL[-2][CECindex + 2 * q.noofCEC][-6] + xdiff,q.statesL[-2][CECindex + 2 * q.noofCEC][-5] + ydiff,q.statesL[-2][CECindex + 2 * q.noofCEC][-4] + zdiff,0,0,0])
              contL.append([q.statesL[-1][CECindex][0],-2,q.statesL[-2][CECindex + q.noofCEC][-6],q.statesL[-2][CECindex + q.noofCEC][-5],q.statesL[-2][CECindex + q.noofCEC][-4],0,0,0])
          else:
              if q.statesL[-1][0] - q.statesL[0][0] < 500000:
                 q.counthopsL[i][1]+=1
                 q.counthopsL[i][3].append(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)
              contL.append([q.statesL[-1][CECindex][0],-2,q.statesL[-2][CECindex + q.noofCEC][-6] + xdiff,q.statesL[-2][CECindex + q.noofCEC][-5] + ydiff,q.statesL[-2][CECindex + q.noofCEC][-4] + zdiff,0,0,0])              
              distL.append([q.statesL[-1][CECindex][0],-3,q.statesL[-2][CECindex + 2 * q.noofCEC][-6],q.statesL[-2][CECindex + 2 * q.noofCEC][-5],q.statesL[-2][CECindex + 2 * q.noofCEC][-4],0,0,0])
        for i in range(0,q.noofCEC):
            q.statesL[-1].append(contL[i])
        for i in range(0,q.noofCEC):
            q.statesL[-1].append(distL[i])

    if q.computeC1Sdist:        
        calcTseContourPlot(q,q.TseL1,q.TseL2,q.TseL3)
    
    count = 0
    for type in q.MSDtypes:
        calcMSD(q,q.bigMSDL[count],type)               
        count+=1
    count = 0
    for type in q.angletypes:
        calcangle(q,q.bigangleL[count],type)               
        count+=1
    if q.calculate1DSqdisp:
        count = 1
        for type in q.Sqdisp:
            calc1DSqdisp(q,q.big1DSqdisp[count],type)               
            count+=1        

    if q.MSDdecomposition:  #reset whether things have reacted or not.
        for item in q.CECorderL:
            item[2] = 0    
                
    q.bigRDFL[0]+=1
    count = 1 #the 0th place is counter
    for pair in q.pairL:                    
        calcRDF(q,q.bigRDFL[count])
        count+=1
    
    count = 0


    if q.calculateclosest == True:
        q.bigClosestL[0]+=1
        q.bigConditionalL[0]+=1
        count = 1  #the 0th place is counter
        for pair in q.closestL:
            calcClosest(q,q.bigClosestL[count], count)
            count+=1
    if q.calculate_endtoenddist:
        q.currentdist=endtoenddist(q,q.currentdist)
    if q.findwaterclusters:
        findwaterclusters(q)
    if q.findFPT:
        findFPT(q)


