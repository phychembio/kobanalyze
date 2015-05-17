#modified on 1/5/2015
#now use default options in kobinfo.py

import math
import sys
import datetime
from koblib import *
from kobspeciallib import *
from kobcompute import *
from kobspecialcompute import *
from kobinfo import*

starttime=datetime.datetime.now()
args=sys.argv[1:]
if len(args)==0:
    args.append("")

q=KobInfo()
q.xl=2


OWID=3
q.OWID=OWID
HWID=4
q.HWID=HWID
OHID=3
q.OHID=OHID
HHID=4
q.HHID=HHID
FID=2
ClID=5
NID=6


q.SID=10
q.C1ID=2
q.OID=5
q.NID=6

prefix="80CsFCl/"

q.trajfile=prefix+'short'+args[0]
q.evbfile=prefix+''+args[0]
q.molfile=prefix+'start.lammpstrj'
driftfile='drift.data'
q.findFPT=False
q.calculate_waterresidence=True

q.MSDtypes=[1]
q.pairL=[]
#q.pairL=[([NID],[ClID]),([NID],[FID]),([NID],[OWID]),([FID],[ClID]),([FID],[OWID]),([ClID],[OWID])]                                                      

if q.findFPT:
    q.FPTbinsize=2000
    q.FPTtypes=[(q.OWID,3.5,4000)] #type,distance threshold, MPFT to decide which are long-stays. Put -999 if it is not known.
    q.FPTpairL=[]
    q.calculateFPTproperties=True    
    if q.calculateFPTproperties:
        q.FPTpairL=[([OWID],[OWID]),([OWID],[FID])]    
    #the way to use this is
    #1) generate the distribution of FPT and find th MFPT
    #2) write out which atoms are long stays (stay > MFPT)
    #3) Using the file that says which are long stays, calculate properties

if q.calculate_residence:
    q.respairL=[([NID],[ClID],6.5,2000),([NID],[FID],6.5,2000),([NID],[FID],8.3,2000)] #([A],[B],Max R_A-B, t*)                                              
#q.respairL=[([NID],[FID],7.25,2000)] #([A],[B],Max R_A-B, t*)  

if q.calculate_waterresidence:
    q.WresOOdist=3.5
    q.WresHOOangle=30
    q.Wreststar=2000
    q.collectblockerdata=True
    if q.collectblockerdata:
        q.blockerL=[[FID,ClID]]

if q.calculateclosest:
    q.closestL=[([7],[q.NID], False)]  #the true/false is for whether to calculate the conditional coordination numbers       
    #q.closestL=[([FID],[NID],True),([ClID],[NID], True),([NID],[NID], False)]  #the true/false is for whether to calculate the conditional coordination numbers
    q.conditionL=[([FID],[OWID],3.25),([ClID],[OWID],3.9)]
    q.notpolymerL=[OWID,HWID,NID,98] #for calculating "closest", these are not included

if q.calculate1DSqdisp:
    q.Sqdisp=[9]


initialize(q)

if q.correctdrift:
    q.driftvL=calcdrift(prefix+driftfile)

stepskipped=0 #time steps already skipped, no need to change


trajf=open(q.trajfile+q.trajext,'rU')         
if q.needCEC:
    evbf=open(q.evbfile+q.evbext,'rU') 
else:
    evbf=None
    
          
    
if q.needCEC:        
    currenttime=eventime(evbf,trajf,q)
    etime=currenttime      
else:
    currenttime=trajtime(trajf,q)  
####End initialize#####    

tempstate=inputstate(trajf,q,currenttime)
specialinit(q)

if q.computeC1Sdist:
    #q.TseL1=q.Typemap[q.SID] 
    #q.TseL2=q.Typemap[q.C1ID]
    #q.TseL3=q.Typemap[q.OID] 
    q.TseL1=q.Typemap[q.NID] #NID   
    
       
    q.TseL2=[]
    indexL=range
    count=0
    noofchains=1
    iterable=iter(q.Typemap[1][1:])
    for index in iterable:                
        if count<20:
            q.TseL2.append(index)
            for i in range(5):
                dummy=next(iterable)
            count+=1            
            
        else:
            if noofchains<4:
                for i in range(1):
                    dummy=next(iterable)
                count=0
                noofchains+=1
            else:
                break         

    q.TseL3=q.Typemap[3][0::6] #CX1ID



if q.needCEC:
        (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,currenttime,q)

if q.calculate_endtoenddist:
    #q.headtype=q.Typemap[1] #PFSA
    #q.tailtype=q.Typemap[4] #PFSA

    ####################PVBTMA for 4 chains
    noofchains=4
    q.headtype=[]
    q.tailtype=[]
    indexmultipler=len(q.Typemap[1])/noofchains   #PVBTMA for 4 chains
    for i in range(noofchains):
        q.headtype.append(q.Typemap[1][i*indexmultipler])
        if i>=1:
           q.tailtype.append(q.Typemap[1][i*indexmultipler-1])
    q.tailtype.append(q.Typemap[1][-1])
    #################### end######PVBTMA for 4 chains
    q.currentdist=[]
    for i in range(len(q.headtype)):
       q.currentdist.append([0,0,0]) #X,X^2,count
                   
fixpbc(tempstate,q)
q.statesL.append(tempstate)

while True:
    try:
        currenttime=trajtime(trajf,q)  
        q.time=currenttime
        tempstate=inputstate(trajf,q,currenttime)  
        if q.needCEC:
            if  q.statesL[-1][0]<currenttime:
                etime=evbtime(evbf)
                (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,etime,q)                          
                #print("%d     %d"%(q.statesL[-1][0],currenttime))
                while etime<currenttime:
                    etime=evbtime(evbf) 
                    (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,etime,q)
                                     
        if q.statesL[-1][0]<currenttime:                             
            fixpbc(tempstate,q)
            q.statesL.append(tempstate)             

####Check/Debug########

        if currenttime%10000==0:
            print currenttime
            sys.stdout.flush()    
        if currenttime%500000==0:
            printresults(args[0],q)
            specialprintresults(args[0],q)


        if q.MSDdecomposition==False and len(q.statesL)>0 and len(q.statesL[-1])!=len(q.statesL[0])\
            or q.MSDdecomposition==True and len(q.statesL)>2 and len(q.statesL[-1])+2*q.noofCEC!=len(q.statesL[0]):
            print "Current Time: " +str(currenttime)
            print "latest state size: " +str(len(q.statesL[-1])) 
            print "initial state size: " +str(len(q.statesL[0]))
            sys.exit("state sizes mismatch..quitting")     
           
####End Check/Debug########


        compute(q) #calculate quantities from statesL        
        specialcompute(q) #calculate for things that are not so common

        uppert=upper(currenttime)
        beforesize=len(q.statesL)

        kobclean(q.statesL,currenttime)  #clean up states

        #if uppert==currenttime:                                                                                           
        
        
    except StopIteration:           
        break


    

#for i in range(len(q.counthopsL)):
#    print("%d %d\n"%(q.counthopsL[i][0],q.counthopsL[i][1]))

printresults(args[0],q)
specialprintresults(args[0],q)

endtime=datetime.datetime.now()
print endtime-starttime


