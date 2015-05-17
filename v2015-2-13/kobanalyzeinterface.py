#modified on 5/23/2013
#added modules koblib and kobcompute

import math
import sys
import datetime
from koblib import *
from kobcompute import *
from kobinfo import*



starttime=datetime.datetime.now()
args=sys.argv[1:]
if len(args)==0:
    args.append("")

q=KobInfo()

OWID=1
q.OWID=OWID
HWID=2
q.HWID=HWID
OHID=3
q.OHID=OHID
HHID=4
q.HHID=HHID
FID=9
NID=6
ClID=10

prefix="interface/"
q.trajfile=prefix+'shortNVT'+args[0]
q.evbfile=prefix+'shortNVT'+args[0]
#trajfile="shortPVBTMA"
#evbfile=""
q.trajext='.lammpstrj'
q.evbext='.evb'
q.molfile=prefix+'start.lammpstrj'
driftfile='drift.data'
q.binsize=0.05
q.ci_binsize=0.025
q.anglebinsize=1.25*(4*math.atan(1))/180
q.skipstep=0

q.needmolID=False 
q.correctdrift=False
#needed for reactive dynamics

q.computeC1Sdist=False
q.MSDdecomposition=False
q.calculateclosest=False
q.calculate1DSqdisp=False
q.c_i_info=False
q.calculate_residence=False
q.calculate_endtoenddist=False


q.interface=True
q.bigionL=[]
q.iontype=[1,3,4,8]
q.list4=[]
q.massL=[16,1.008,35.45,12.01,1.008,12.01,1.008,14.01]
q.initialzcenterofmass=-1


q.xl=2
q.SID=10
q.C1ID=2
q.OID=5


q.MSDtypes=[]
q.angletypes=[]

q.pairL=[]




if q.calculate_residence:
    q.respairL=[([NID],[ClID],6.5,2000),([NID],[FID],6.5,2000),([NID],[FID],8.3,2000)] #([A],[B],Max R_A-B, t*)
#q.respairL=[([NID],[FID],7.25,2000)] #([A],[B],Max R_A-B, t*)

if q.calculateclosest:
    q.closestL=[([FID],[NID],True),([ClID],[NID], True),([NID],[NID], False)]  #the true/false is for whether to calculate the conditional coordination numbers
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


if q.calculate_residence:
    setupbigtrackingresL(q)
if q.needCEC:
        (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,currenttime,q)
                   
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


        if q.MSDdecomposition==False and len(q.statesL)>0 and len(q.statesL[-1])!=len(q.statesL[0])\
            or q.MSDdecomposition==True and len(q.statesL)>2 and len(q.statesL[-1])+2*q.noofCEC!=len(q.statesL[0]):
            print "Current Time: " +str(currenttime)
            print "latest state size: " +str(len(q.statesL[-1])) 
            print "initial state size: " +str(len(q.statesL[0]))
            sys.exit("state sizes mismatch..quitting")     
           
####End Check/Debug########


        interfacecompute(q) #calculate quantities from statesL

        uppert=upper(currenttime)
        beforesize=len(q.statesL)

        kobclean(q.statesL,currenttime)  #clean up states

        #if uppert==currenttime:                                                                                           
        
        
    except StopIteration:           
        break


    

#for i in range(len(q.counthopsL)):
#    print("%d %d\n"%(q.counthopsL[i][0],q.counthopsL[i][1]))

printresults(args[0],q)


endtime=datetime.datetime.now()
print endtime-starttime


