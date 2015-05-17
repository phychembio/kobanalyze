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


def calselfsqdist(statesL,previndex,currindex,typeid,q):
    Typemap=q.Typemap
    noofatom=statesL[currindex][1]    
    sqsum=0
    
    timediff=statesL[currindex][0]-statesL[previndex][0]
    if q.correctdrift==False:
        for i in Typemap[typeid]:
              xdiff=statesL[currindex][i][2]-statesL[previndex][i][2]+(statesL[currindex][i][5]-statesL[previndex][i][5])*q.xboxwidth
              ydiff=statesL[currindex][i][3]-statesL[previndex][i][3]+(statesL[currindex][i][6]-statesL[previndex][i][6])*q.yboxwidth
              zdiff=statesL[currindex][i][4]-statesL[previndex][i][4]+(statesL[currindex][i][7]-statesL[previndex][i][7])*q.zboxwidth
              sqsum+=xdiff**2+ydiff**2+zdiff**2    
    else:
        for i in Typemap[typeid]:
              xdiff=statesL[currindex][i][2]-statesL[previndex][i][2]+(statesL[currindex][i][5]-statesL[previndex][i][5])*q.xboxwidth-q.driftvL[0]*timediff
              ydiff=statesL[currindex][i][3]-statesL[previndex][i][3]+(statesL[currindex][i][6]-statesL[previndex][i][6])*q.yboxwidth-q.driftvL[1]*timediff
              zdiff=statesL[currindex][i][4]-statesL[previndex][i][4]+(statesL[currindex][i][7]-statesL[previndex][i][7])*q.zboxwidth-q.driftvL[2]*timediff
              sqsum+=xdiff**2+ydiff**2+zdiff**2
    
    if sqsum >len(Typemap[typeid])*6*0.005*timediff:
          print("current: %d previous: %d  , No. of Atoms of this type: %d , sqsum: %f"%(statesL[currindex][0],statesL[previndex][0],len(Typemap[typeid]),sqsum))                 
              
    return (sqsum,len(Typemap[typeid]))

def findwhichtimeslot(timediff,smallesttime):
    if timediff==0:
        return -1
    else:
        time_unit=timediff/smallesttime #e.g timediff=400000, smallesttime =1000, time_unit=400
        whichpow=math.floor(math.log10(time_unit)) #whichpow=2
        firstdigit=time_unit/10**whichpow #firstdigit=4
        return int(whichpow*10+firstdigit-(whichpow+1)) #return 23

def calclastMSD(q,type):
    sum=0
    sqdistinfo=calselfsqdist(q.statesL,0,-1,type,q)
    timediff=q.statesL[-1][0]-q.statesL[0][0]
    return (timediff,sqdistinfo[0]*1./sqdistinfo[1])


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
              sqdistinfo=calselfsqdist(q.statesL,i,-1,type,q)
              if timediff==currtimediff:                  
                  MSDL.append([timediff,sqdistinfo[0],sqdistinfo[1]]) #sqdistinfo[0]=sqdist, [1]=count
              elif timediff<currtimediff:
                  slot=findwhichtimeslot(timediff,q.smalleststep)                
                  MSDL[slot+1][1]+=sqdistinfo[0] #the +1 is for the label
                  MSDL[slot+1][2]+=sqdistinfo[1]                        
     
                  
                               
def calcResT(q,trackingresL,resL):
       pair=resL[0]       
       list1=q.Typemap[pair[0][0]]
       list2=q.Typemap[pair[1][0]]
       maxdist=pair[2]
       t_star=pair[3]
       currentt=q.statesL[-1][0]
       for item in pair[0][1:]:
          list1+=q.Typemap[item]
       for item in pair[1][1:]:
           list2+=q.Typemap[item]

       pos1=0                           
       for i in list1:            
            pos2=0
            for j in list2:              
                if j!=i:                
                    xdiff1=minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],q.xboxwidth)
                    ydiff1=minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],q.yboxwidth)
                    zdiff1=minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],q.zboxwidth)
                    dist1=math.sqrt(xdiff1**2+ydiff1**2+zdiff1**2)                                        
                    #if q.statesL[-1][i][0]==1341 and q.statesL[-1][j][0]==5878:                  
                    #    if dist1<maxdist:
                    #        q.debugrestf2.write(" %d  , %d , %f , %d, YES\n "%(q.statesL[-1][i][0],q.statesL[-1][j][0],dist1, currentt,))
                    #    else:
                    #        q.debugrestf2.write(" %d  , %d , %f , %d, NO\n "%(q.statesL[-1][i][0],q.statesL[-1][j][0],dist1, currentt,))
                startt=trackingresL[pos1][pos2][0]   
                lastt=trackingresL[pos1][pos2][1]        
                if dist1<maxdist:                           
                    if startt==-1:                         
                         trackingresL[pos1][pos2][0]=currentt
                         trackingresL[pos1][pos2][1]=currentt
                    else:                         
                         if currentt-(lastt+q.smalleststep)< t_star: # note: dist1<maxdist already confirms there is residence in this step.
                                                                     # see the difference about 20 lines later
                            trackingresL[pos1][pos2][1]=currentt
                         else: #new code doesn't have this part. I am stil confused (2/20/2015) why it existed.
#================================================================                  
                            if startt>q.statesL[0][0]+t_star: # remember that the first time this script actually reads is the second time slot.
                                timediff=(lastt+ q.smalleststep)-startt                 
                                maxbin=len(resL)-1 
                                whichbin=timediff/q.smalleststep
                                if whichbin>maxbin:                      
                                    while maxbin<whichbin: # -1 is for the label for the pair at the first slot                          
                                        if len(resL)==1:
                                            resL.append([timediff,0])
                                        else:
                                            resL.append([resL[-1][0]+q.smalleststep,0])  
                                        maxbin+=1
                                    resL[-1][1]+=1 #last entry
                                else:                       
                                    resL[whichbin][1]+=1   
                                #q.debugrestf.write("%d, %d , %d, %d , %d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep))
#================================================================                        
                            trackingresL[pos1][pos2][0]=currentt
                            trackingresL[pos1][pos2][1]=currentt
                else:                      
                    if lastt!=-1 and currentt-lastt >= t_star:        
                        # note: dist1>=maxdist already confirms there is no residence in this time.
#================================================================
                        if startt>q.statesL[0][0]+t_star: # remember that the first time this script actually reads is the second time slot.
                            timediff=(lastt+ q.smalleststep)-startt                                
                            maxbin=len(resL)-1 
                            whichbin=timediff/q.smalleststep
                            if whichbin>maxbin:                      
                                while maxbin<whichbin: # -1 is for the label for the pair at the first slot                          
                                    if len(resL)==1:
                                        resL.append([timediff,0])
                                    else:
                                        resL.append([resL[-1][0]+q.smalleststep,0])  
                                    maxbin+=1
                                resL[-1][1]+=1 #last entry
                            else:                       
                                resL[whichbin][1]+=1   
                            #q.debugrestf.write("%d, %d , %d, %d , %d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep))
#================================================================                
                        trackingresL[pos1][pos2][0]=-1
                        trackingresL[pos1][pos2][1]=-1
                pos2+=1
            pos1+=1
       
def construct_resT_CF(q,resL):
    L=[]
    for i in range(len(resL)):
        L.append([i*q.smalleststep,0])
    L[0][1]=resL[0][1]    
    for i in range(1,len(resL)):
        for j in range(i,-1,-1):
            L[j][1]+=resL[i][1]
    N=L[0][1]
    for item in L:
        item[1]=item[1]*1./N
    return L

def calc_meanresT(resL):
    Normalization=0
    sum=0
    for item in resL:
        sum+=item[0]*item[1]
        Normalization+=item[1]
    return sum*1./Normalization
     
                                         

def calcRDF(q,RDFL):
    pair=RDFL[0]    
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
    
                xdiff=minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],q.xboxwidth)
                ydiff=minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],q.yboxwidth)
                zdiff=minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],q.zboxwidth)
                dist=math.sqrt(xdiff**2+ydiff**2+zdiff**2)
                

                #  print("%f  %f  %f\n"%(q.statesL[-1][j][2],q.statesL[-1][i][2],xdiff))
                #  print("%f  %f  %f\n"%(q.statesL[-1][j][3],q.statesL[-1][i][3],ydiff))
                #  print("%f  %f  %f\n"%(q.statesL[-1][j][4],q.statesL[-1][i][4],zdiff))
                #  print ("%f \n\n"%dist)
                
                whichbin=int(dist/q.binsize)+1 # +1 to account for the label slot            
                koblib.histo(RDFL,q.binsize,whichbin,1)
   

def calcClosest(q,ClosestL,listcount):
    pair=ClosestL[0]
    list1=q.Typemap[pair[0][0]]    
    for item in pair[0][1:]:
        list1+=q.Typemap[item]

    list2=[]
    if pair[1][0]!=999:
        list2=q.Typemap[pair[1][0]]
        for item in pair[1][1:]:
            list2+=q.Typemap[item]
    else:        
        for i in range(1,len(q.Typemap)-1):
            if i not in notpolymerL:
                list2=list2+q.Typemap[i]         
        if needCEC==True: #if CEC, the last atom type is CEC and not polymer
            list2=list2+q.Typemap[-1]         
        elif (i+1) not in notpolymerL:
                list2=list2+q.Typemap[-1]            
            
        
    
    closestindex=-1
    for i in list1:            
            mindist=q.zboxwidth            
            for j in list2:
                 if j!=i:                        
                       xdiff=minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],q.xboxwidth)
                       ydiff=minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],q.yboxwidth)
                       zdiff=minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],q.zboxwidth)
                       dist=math.sqrt(xdiff**2+ydiff**2+zdiff**2)
                       if dist<mindist:
                              mindist=dist                                                 
            whichbin=int(mindist/q.binsize)+1 # +1 to account for the label slot
            koblib.histo(ClosestL,q.binsize,whichbin,1)
            if pair[2]==True:                
                list3=q.Typemap[q.conditionL[listcount-1][1][0]]
                cutoff=q.conditionL[listcount-1][2]
                atomcount=0
                for k in list3:
                    if k!=i:
                       xdiff=minimagedist(q.statesL[-1][i][2],q.statesL[-1][k][2],q.xboxwidth)
                       ydiff=minimagedist(q.statesL[-1][i][3],q.statesL[-1][k][3],q.yboxwidth)
                       zdiff=minimagedist(q.statesL[-1][i][4],q.statesL[-1][k][4],q.zboxwidth)
                       dist=math.sqrt(xdiff**2+ydiff**2+zdiff**2)
                       if dist<cutoff:
                           atomcount+=1

                binL=q.bigConditionalL[listcount]
                maxbin=len(binL)-1 
                if whichbin>maxbin:                      
                        while maxbin<whichbin: # -1 is for the label for the pair at the first slot                          
                            if len(binL)==1:
                              binL.append([q.binsize/2,0,0])
                            else:
                              binL.append([binL[-1][0]+q.binsize,0,0])  
                            maxbin+=1
                        binL[-1][1]+=atomcount #last entry
                        binL[-1][2]+=1 #last entry
                else:
                      binL[whichbin][1]+=atomcount    
                      binL[whichbin][2]+=1

        

count1=0
count2=0

def calcC1Sclosest(q):
            
    count=0
    for i in q.Typemap[q.SID]:
        
        j=q.Typemap[q.C1ID][count]
        k=q.Typemap[q.OID][count]
        
        
        xdiff1=minimagedist(q.statesL[-1][i][2],q.statesL[-1][j][2],q.xboxwidth)
        ydiff1=minimagedist(q.statesL[-1][i][3],q.statesL[-1][j][3],q.yboxwidth)
        zdiff1=minimagedist(q.statesL[-1][i][4],q.statesL[-1][j][4],q.zboxwidth)
        dist=math.sqrt(xdiff1**2+ydiff1**2+zdiff1**2)
        
        xdiff2=minimagedist(q.statesL[-1][k][2],q.statesL[-1][j][2],q.xboxwidth)
        ydiff2=minimagedist(q.statesL[-1][k][3],q.statesL[-1][j][3],q.yboxwidth)
        zdiff2=minimagedist(q.statesL[-1][k][4],q.statesL[-1][j][4],q.zboxwidth)
        
        temp1=[xdiff1,ydiff1,zdiff1]
        len1=koblib.lenvec(temp1)
        SC1vec=[xdiff1/len1,ydiff1/len1,zdiff1/len1]
        temp2=[xdiff2,ydiff2,zdiff2]
        len2=koblib.lenvec(temp2)
        Cb1C2vec=[xdiff2/len2,ydiff2/len2,zdiff2/len2]
       
        whichbin=int(dist/q.binsize)+1 # +1 to account for the label slot
        
        angle=math.acos(koblib.dotvec(SC1vec,Cb1C2vec))
        whichanglebin=int(angle/q.anglebinsize)
        maxbin=len(q.C1SL)-1
        if whichanglebin <80:
          if whichbin>maxbin:
            while maxbin<whichbin: # -1 is for the label for the pair at the first slot                          
                if len(q.C1SL)==1:
                    q.C1SL.append([q.binsize/2,0])
                      
                else:
                    q.C1SL.append([q.C1SL[-1][0]+q.binsize,0])                     
                    maxbin+=1
                for i in range(0,80):
                    q.C1SL[-1].append(0)
            q.C1SL[-1][1]+=1 #last entry
            q.C1SL[-1][1+1+whichanglebin]+=1 #last entry
          else:
            q.C1SL[whichbin][1]+=1
            q.C1SL[whichbin][1+1+whichanglebin]+=1
        else:
           print whichanglebin
           print "Not collected"
        count+=1

def collect_ci(ciL,q,noofstates):
    q.bigciL[0]+=1
    for i in range(2):              
        cisq=ciL[i]**2
        whichbin=int(cisq/q.ci_binsize)+1 # +1 to account for the label slot
        koblib.histo(q.bigciL[i+1],q.ci_binsize,whichbin,1)         
    ci_sqdiff=ciL[0]*ciL[0]-ciL[1]*ciL[1]
    koblib.histo(q.bigciL[3],q.ci_binsize,whichbin,-math.log(ciL[0]/ciL[1]))  
    whichbin=int(ci_sqdiff/q.ci_binsize)+1
    whichbin=noofstates+1    
    koblib.histo(q.NstateL,1,whichbin,1)
    

def compute(q):         
    if q.MSDdecomposition:
        noofatom=q.statesL[-1][1]        
        distL=[]
        contL=[]
        if len(q.statesL[0])==3+noofatom+q.noofCEC:
            for i in range(0,q.noofCEC):
                print i
                CECindex=noofatom+3+i
                q.statesL[0].append([q.statesL[-1][CECindex][0],-2,0,0,0,0,0,0])
            for i in range(0,q.noofCEC):
                CECindex=noofatom+3+i
                q.statesL[0].append([q.statesL[-1][CECindex][0],-3,0,0,0,0,0,0])
        for i in range(0,q.noofCEC):
          CECindex=noofatom+3+i
          xdiff=q.statesL[-1][CECindex][2]-q.statesL[-2][CECindex][2]+(q.statesL[-1][CECindex][5]-q.statesL[-2][CECindex][5])*q.xboxwidth
          ydiff=q.statesL[-1][CECindex][3]-q.statesL[-2][CECindex][3]+(q.statesL[-1][CECindex][6]-q.statesL[-2][CECindex][6])*q.yboxwidth
          zdiff=q.statesL[-1][CECindex][4]-q.statesL[-2][CECindex][4]+(q.statesL[-1][CECindex][7]-q.statesL[-2][CECindex][7])*q.zboxwidth

          reacted=False
          CEC2ndorderL=sorted(q.CECorderL, key=lambda k:k[1])
          #if q.CECorderL[i][2]==1:
          if CEC2ndorderL[i][2]==1:
              reacted=True
          q.hopsL[i].append(q.statesL[-1][CECindex][0])
          if reacted==True:
              if q.statesL[-1][0]-q.statesL[0][0]<500000:
                 q.counthopsL[i][0]+=1
                 q.counthopsL[i][2].append(xdiff**2+ydiff**2+zdiff**2)
              distL.append([q.statesL[-1][CECindex][0],-3,q.statesL[-2][CECindex+2*q.noofCEC][-6]+xdiff,q.statesL[-2][CECindex+2*q.noofCEC][-5]+ydiff,q.statesL[-2][CECindex+2*q.noofCEC][-4]+zdiff,0,0,0])
              contL.append([q.statesL[-1][CECindex][0],-2,q.statesL[-2][CECindex+q.noofCEC][-6],q.statesL[-2][CECindex+q.noofCEC][-5],q.statesL[-2][CECindex+q.noofCEC][-4],0,0,0])
          else:
              if q.statesL[-1][0]-q.statesL[0][0]<500000:
                 q.counthopsL[i][1]+=1
                 q.counthopsL[i][3].append(xdiff**2+ydiff**2+zdiff**2)
              contL.append([q.statesL[-1][CECindex][0],-2,q.statesL[-2][CECindex+q.noofCEC][-6]+xdiff,q.statesL[-2][CECindex+q.noofCEC][-5]+ydiff,q.statesL[-2][CECindex+q.noofCEC][-4]+zdiff,0,0,0])              
              distL.append([q.statesL[-1][CECindex][0],-3,q.statesL[-2][CECindex+2*q.noofCEC][-6],q.statesL[-2][CECindex+2*q.noofCEC][-5],q.statesL[-2][CECindex+2*q.noofCEC][-4],0,0,0])
        for i in range(0,q.noofCEC):
            q.statesL[-1].append(contL[i])
        for i in range(0,q.noofCEC):
            q.statesL[-1].append(distL[i])

    if q.computeC1Sdist:        
        calcC1Sclosest(q)
    
    count=0
    for type in q.MSDtypes:
        calcMSD(q,q.bigMSDL[count],type)               
        count+=1

    if q.MSDdecomposition:  #reset whether things have reacted or not.
        for item in q.CECorderL:
            item[2]=0    
                
    q.bigRDFL[0]+=1
    count=1 #the 0th place is counter
    for pair in q.pairL:                    
        calcRDF(q,q.bigRDFL[count])
        count+=1
    
    count=0
    for pair in q.respairL:
        calcResT(q,q.bigtrackingresL[count],q.bigresL[count])
        count+=1

    if q.calculateclosest==True:
        q.bigClosestL[0]+=1
        q.bigConditionalL[0]+=1
        count=1  #the 0th place is counter
        for pair in q.closestL:
            calcClosest(q,q.bigClosestL[count], count)
            count+=1
  