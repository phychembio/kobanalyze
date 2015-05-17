#!/usr/bin/python
import re
import sys
import math
import koblib

def upper(s):
    multipleofsmallest=s
    done = False
    counter = 0
    temp=multipleofsmallest
    while not done and s!=0:
        if temp%10 !=0:
            done = True
        else:
            temp=temp/10
            counter+=1 
    return multipleofsmallest- temp%10*(10**counter) + 10**(counter+1)

def kobclean(statesL,currenttime):
    for state in statesL[1:]:
            s=state[0]
            if upper(s)==currenttime:
                statesL.remove(state)                                     

        

def out_C_sq(filename, starttime):
    f=open(filename,'rU')
    timeskipto(f,starttime)
    fout=open('c_sq.txt','w')
 
    timestep=starttime

    for line in f:
        match0=re.match("TIMESTEP (\d+)",line)
        match=re.match("EIGEN_VECTOR",line)
        if match0:
            timestep=match0.group(1)
 
        if match and int(timestep)%100==0:
           L=[]
           for c in f.next().split():
             L.append(float(c)*float(c))
           L=sorted(L, reverse=True)
           fout.write("%f   %f \n"% (L[0], L[1]))
           
           
def minimagedist(r2,r1,boxwidth):
    diff=r2-r1
    if math.fabs(diff) > boxwidth/2:
        if diff>0:
            return diff-boxwidth
        else:
            return diff+boxwidth
    else:
        return diff


def calselfsqdist(statesL,previndex,currindex,typeid,Typemap):
   
    noofatom=statesL[currindex][1]
    boxwidth=statesL[currindex][2]*2
    sqsum=0
    count=0    
    
    for i in Typemap[typeid]:
          xdiff=statesL[currindex][i][2]-statesL[previndex][i][2]+(statesL[currindex][i][5]-statesL[previndex][i][5])*boxwidth
          ydiff=statesL[currindex][i][3]-statesL[previndex][i][3]+(statesL[currindex][i][6]-statesL[previndex][i][6])*boxwidth
          zdiff=statesL[currindex][i][4]-statesL[previndex][i][4]+(statesL[currindex][i][7]-statesL[previndex][i][7])*boxwidth
          sqsum+=xdiff**2+ydiff**2+zdiff**2
          count+=1
        
          if sqsum >6500:
              print("current: %d previous: %d  index: %d  sqsum: %f"%(statesL[currindex][0],statesL[previndex][0],i,sqsum))                 
              
    return (sqsum,count)

def findwhichtimeslot(timediff,smallesttime):
    if timediff==0:
        return -1
    else:
        time_unit=timediff/smallesttime #e.g timediff=400000, smallesttime =1000, time_unit=400
        whichpow=math.floor(math.log10(time_unit)) #whichpow=2
        firstdigit=time_unit/10**whichpow #firstdigit=4
        return int(whichpow*10+firstdigit-(whichpow+1)) #return 23


def calcMSD(q,MSDL,type):   
    no_of_timeslots=len(q.statesL)    
    currtimediff=q.statesL[-1][0]-q.statesL[0][0]    
    if no_of_timeslots>0:
         for i in range(no_of_timeslots-2,-1,-1): #reading the times backards
            timediff=q.statesL[-1][0]-q.statesL[i][0]
            print("%d:%d %d:%d   %d"%(no_of_timeslots-1,q.statesL[-1][0],i, q.statesL[i][0],timediff))            
                    
            if timediff>0:
                basetime=int(10**math.floor(math.log10(timediff)))
            else:
                sys.exit("Error: time difference >= 0") 
                

            if timediff%basetime==0:           
              sqdistinfo=calselfsqdist(q.statesL,i,-1,type,q.Typemap)
              if timediff==currtimediff:                  
                  MSDL.append([timediff,sqdistinfo[0],sqdistinfo[1]]) #sqdistinfo[0]=sqdist, [1]=count
              elif timediff<currtimediff:
                  slot=findwhichtimeslot(timediff,q.smalleststep)                
                  MSDL[slot+1][1]+=sqdistinfo[0] #the +1 is for the label
                  MSDL[slot+1][2]+=sqdistinfo[1]                        
                  
def calcResT(q,resL,pair,maxdist):

    no_of_timeslots=len(q.statesL)        
    currtimediff=q.statesL[-1][0]-q.statesL[0][0]    
    if no_of_timeslots>0:
         for t_0 in range(no_of_timeslots-1,-1,-1): #reading the times backards, the second index=-1 means it starts with the last element, the their index=-1 means it stops at 0!
            timediff=q.statesL[-1][0]-q.statesL[t_0][0]
            if timediff>0:
                basetime=int(10**math.floor(math.log10(timediff)))
            
            if timediff==0 or timediff%basetime==0:                   
                boxwidth=q.statesL[-1][2]*2
                list1=q.Typemap[pair[0][0]]
                list2=q.Typemap[pair[1][0]]
                for item in pair[0][1:]:
                    list1+=q.Typemap[item]
                for item in pair[1][1:]:
                    list2+=q.Typemap[item]
                    
                for i in list1:
                        #print i            
                        for j in list2:              
                          if j!=i:                             
                            xdiff1=minimagedist(q.statesL[-1][j][2],q.statesL[t_0][i][2],boxwidth)
                            ydiff1=minimagedist(q.statesL[-1][j][3],q.statesL[t_0][i][3],boxwidth)
                            zdiff1=minimagedist(q.statesL[-1][j][4],q.statesL[t_0][i][4],boxwidth)
                            dist1=math.sqrt(xdiff1**2+ydiff1**2+zdiff1**2)            
                            
                            if dist1<maxdist:
                                xdiff2=minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],boxwidth)
                                ydiff2=minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],boxwidth)
                                zdiff2=minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],boxwidth)
                                dist2=math.sqrt(xdiff2**2+ydiff2**2+zdiff2**2)            

                                if dist2<maxdist:
                                     latest_t_inL=resL[-1][0]
                                     if timediff>latest_t_inL: 
                                         resL.append([timediff,1])
                                     elif timediff<=currtimediff:                                
                                            slot=findwhichtimeslot(timediff,q.smalleststep)  
                                            resL[slot+2][1]+=1 #the +2 in "slot+1" is for the label and the timediff=0 slot
                 
                                 
                                         

def calcRDF(q,RDFL,pair):
    
    boxwidth=q.statesL[-1][2]*2
    list1=q.Typemap[pair[0][0]]
    list2=q.Typemap[pair[1][0]]
    for item in pair[0][1:]:
        list1+=q.Typemap[item]
    for item in pair[1][1:]:
        list2+=q.Typemap[item]


    for i in list1:
            #print i            
            for j in list2:              
              if j!=i:                         
    
                xdiff=minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],boxwidth)
                ydiff=minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],boxwidth)
                zdiff=minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],boxwidth)
                dist=math.sqrt(xdiff**2+ydiff**2+zdiff**2)
                

                #  print("%f  %f  %f\n"%(q.statesL[-1][j][2],q.statesL[-1][i][2],xdiff))
                #  print("%f  %f  %f\n"%(q.statesL[-1][j][3],q.statesL[-1][i][3],ydiff))
                #  print("%f  %f  %f\n"%(q.statesL[-1][j][4],q.statesL[-1][i][4],zdiff))
                #  print ("%f \n\n"%dist)
                
                whichbin=int(dist/q.binsize)+1 # +1 to account for the label slot            
                koblib.histo(RDFL,q.binsize,whichbin)
   

def calcClosest(statesL,pair, notpolymerL, distanceL,Typemap,binsize):
        
    list1=Typemap[pair[0][0]]    
    for item in pair[0][1:]:
        list1+=Typemap[item]

    list2=[]
    if pair[1][0]!=999:
        list2=Typemap[pair[1][0]]
        for item in pair[1][1:]:
            list2+=Typemap[item]
    else:        
        for i in range(1,len(Typemap)-1):
            if i not in notpolymerL:
                list2=list2+Typemap[i]         
        if needCEC==True: #if CEC, the last atom type is CEC and not polymer
            list2=list2+Typemap[-1]         
        elif (i+1) not in notpolymerL:
                list2=list2+Typemap[-1]            
            
        
    boxwidth=statesL[-1][2]*2    
    
    for i in list1:            
            mindist=boxwidth            
            for j in list2:
                 if j!=i:                        
                       xdiff=minimagedist(statesL[-1][j][2],statesL[-1][i][2],boxwidth)
                       ydiff=minimagedist(statesL[-1][j][3],statesL[-1][i][3],boxwidth)
                       zdiff=minimagedist(statesL[-1][j][4],statesL[-1][i][4],boxwidth)
                       dist=math.sqrt(xdiff**2+ydiff**2+zdiff**2)
                       if dist<mindist:
                              mindist=dist         
            
            whichbin=int(dist/q.binsize)+1 # +1 to account for the label slot
            koblib.histo(distanceL,q.binsize,whichbin)
              
 

count1=0
count2=0

def calcC1Sclosest(statesL,pair, C1SL,binsize,atommapC1S):
    maxbin=len(C1SL)-1    
    noofatom=statesL[-1][1]
    boxwidth=statesL[-1][2]*2
    C1SL[0]+=1
    orderedstate=sorted(statesL[-1][3:noofatom+3],key=lambda k: k[0])
    
    if len(atommapC1S)==0:        
        Cb1index=None        
        O1index=None
        Sindex=None
        count=0
        for item in orderedstate:
            if  item[1]==C1ID:
                C1index=count
            elif item[1]==O1ID:
                O1index=count
            elif item[1]==SID:
                atommapC1S.append([count,C1index,O1index])
            count+=1
            

    for item in atommapC1S:
        i=item[0] #S
        j=item[1] #C1
        
        xdiff1=minimagedist(orderedstate[i][2],orderedstate[j][2],boxwidth)
        ydiff1=minimagedist(orderedstate[i][3],orderedstate[j][3],boxwidth)
        zdiff1=minimagedist(orderedstate[i][4],orderedstate[j][4],boxwidth)
        dist=math.sqrt(xdiff1**2+ydiff1**2+zdiff1**2)

        k=item[2] #O1        
        xdiff2=minimagedist(orderedstate[k][2],orderedstate[j][2],boxwidth)
        ydiff2=minimagedist(orderedstate[k][3],orderedstate[j][3],boxwidth)
        zdiff2=minimagedist(orderedstate[k][4],orderedstate[j][4],boxwidth)
        
        temp1=[xdiff1,ydiff1,zdiff1]
        len1=lenvec(temp1)
        SC1vec=[xdiff1/len1,ydiff1/len1,zdiff1/len1]
        temp2=[xdiff2,ydiff2,zdiff2]
        len2=lenvec(temp2)
        Cb1C2vec=[xdiff2/len2,ydiff2/len2,zdiff2/len2]
       
        whichbin=int(dist/binsize)+1 # +1 to account for the label slot
        global count1,count2
        if dotvec(SC1vec,Cb1C2vec)>=0:
            count1+=1
        else:
            count2+=1
        angle=math.acos(dotvec(SC1vec,Cb1C2vec))
        whichanglebin=int(angle/(2.2222222*math.atan(1)/80))
        if whichanglebin <80:
          if whichbin>maxbin:
            while maxbin<whichbin: # -1 is for the label for the pair at the first slot                          
                if len(C1SL)==1:
                    C1SL.append([binsize/2,0])
                      
                else:
                    C1SL.append([C1SL[-1][0]+binsize,0])                     
                    maxbin+=1
                for i in range(0,80):
                    C1SL[-1].append(0)
            C1SL[-1][1]+=1 #last entry
            C1SL[-1][1+1+whichanglebin]+=1 #last entry
          else:
            C1SL[whichbin][1]+=1
            C1SL[whichbin][1+1+whichanglebin]+=1
        else:
           print whichanglebin
           print "Not collected"

def collect_ci(ciL,q,noofstates):
    q.bigciL[0]+=1
    for i in range(1,3):
        ci_binsize=0.05                
        whichbin=int(ciL[i-1]/ci_binsize)+1 # +1 to account for the label slot
        koblib.histo(ciL[i-i],ci_binsize,whichbin)         
    
    whichbin=noofstates+1    
    koblib.histo(q.NstateL,binsize,whichbin)
    

def compute(q):
         
    if q.MSDdecomposition:
        noofatom=q.statesL[-1][1]
        boxwidth=q.statesL[-1][2]*2
        distL=[]
        contL=[]
        if len(q.statesL[0])==3+noofatom+noofCEC:
            for i in range(0,noofCEC):
                print i
                CECindex=noofatom+3+i
                q.statesL[0].append([q.statesL[-1][CECindex][0],-2,0,0,0,0,0,0])
            for i in range(0,noofCEC):
                CECindex=noofatom+3+i
                q.statesL[0].append([q.statesL[-1][CECindex][0],-3,0,0,0,0,0,0])
        for i in range(0,noofCEC):
          CECindex=noofatom+3+i
          xdiff=q.statesL[-1][CECindex][2]-q.statesL[-2][CECindex][2]+(q.statesL[-1][CECindex][5]-q.statesL[-2][CECindex][5])*boxwidth
          ydiff=q.statesL[-1][CECindex][3]-q.statesL[-2][CECindex][3]+(q.statesL[-1][CECindex][6]-q.statesL[-2][CECindex][6])*boxwidth
          zdiff=q.statesL[-1][CECindex][4]-q.statesL[-2][CECindex][4]+(q.statesL[-1][CECindex][7]-q.statesL[-2][CECindex][7])*boxwidth

          reacted=False
          if q.lastorderL[i][2]==1:
              reacted=True
          hopsL[i].append(q.statesL[-1][CECindex][0])
          if reacted==True:
              if q.statesL[-1][0]-q.statesL[0][0]<500000:
                 counthopsL[i][0]+=1
                 counthopsL[i][2].append(xdiff**2+ydiff**2+zdiff**2)
              distL.append([q.statesL[-1][CECindex][0],-3,q.statesL[-2][CECindex+2*noofCEC][-6]+xdiff,q.statesL[-2][CECindex+2*noofCEC][-5]+ydiff,q.statesL[-2][CECindex+2*noofCEC][-4]+zdiff,0,0,0])
              contL.append([q.statesL[-1][CECindex][0],-2,q.statesL[-2][CECindex+noofCEC][-6],q.statesL[-2][CECindex+noofCEC][-5],q.statesL[-2][CECindex+noofCEC][-4],0,0,0])
          else:
              if q.statesL[-1][0]-q.statesL[0][0]<500000:
                 counthopsL[i][1]+=1
                 counthopsL[i][3].append(xdiff**2+ydiff**2+zdiff**2)
              contL.append([q.statesL[-1][CECindex][0],-2,q.statesL[-2][CECindex+noofCEC][-6]+xdiff,q.statesL[-2][CECindex+noofCEC][-5]+ydiff,q.statesL[-2][CECindex+noofCEC][-4]+zdiff,0,0,0])              
              distL.append([q.statesL[-1][CECindex][0],-3,q.statesL[-2][CECindex+2*noofCEC][-6],q.statesL[-2][CECindex+2*noofCEC][-5],q.statesL[-2][CECindex+2*noofCEC][-4],0,0,0])
        for i in range(0,noofCEC):
            q.statesL[-1].append(contL[i])
        for i in range(0,noofCEC):
            q.statesL[-1].append(distL[i])

    if q.computeC1Sdist:
        calcC1Sclosest(q)
    
    count=0
    for type in q.MSDtypes:
        calcMSD(q,q.bigMSDL[count],type)               
        count+=1
               
    q.bigRDFL[0]+=1
    count=1 #the 0th place is counter
    for pair in q.pairL:                    
        calcRDF(q,q.bigRDFL[count],pair)
        count+=1
    
    count=0
    for pair in q.respairL:
        calcResT(q,q.bigresL[count],pair,pair[2])

    if q.calculateclosest==True:
        q.bigClosestL[0]+=1
        count=1  #the 0th place is counter
        for pair in closestL:
            calcClosest(q,q.bigClosestL[count],pair)
            count+=1
    