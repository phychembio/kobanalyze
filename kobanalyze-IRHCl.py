#!/usr/bin/env python

#modified on 1/16/2014
#vibrational spectrum possible
#added modules koblib and kobcompute

import re
import math
import sys

from koblib import *
from kobcompute import *
from kobinfo import*
from operator import itemgetter

args=sys.argv[1:]
if len(args)==0:
    args.append("")




prefix="IR/"
q=KobInfo()

q.OWID=OWID=1
q.HWID=HWID=2
q.OHID=OHID=3
q.HHID=HHID=4
q.ClID=ClID=5
#q.partialchargeL=[99999, -0.82,0.41]
q.partialchargeL=[99999,-0.82,0.41,-0.5,0.5,-1]
q.trajfile=prefix+'shortIR'+args[0]
q.evbfile=prefix+'shortIR'+args[0]
#trajfile="shortPVBTMA"
#evbfile=""
q.trajext='.lammpstrj'
q.evbext='.evb'


q.IR=True #whether to calculate IR spectrum



driftfile='drift.data'
q.binsize=0.05
q.ci_binsize=0.025
q.anglebinsize=1.25*(4*math.atan(1))/180
q.skipstep=0

q.molfile=prefix+'start.lammpstrj' 
q.needmolID=False #needed for reactive dynamics
q.correctdrift=False


q.computeC1Sdist=False
q.MSDdecomposition=False
q.calculateclosest=False
q.calculate1DSqdisp=False
q.c_i_info=False
q.calculate_residence=False
q.calculate_endtoenddist=False


q.MSDtypes=[-1]
q.pairL=[]
#q.pairL=[([NID],[ClID]),([NID],[FID]),([NID],[OWID]),([FID],[ClID]),([FID],[OWID]),([ClID],[OWID])]                                                      



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
setupTypemap(tempstate,q,q.Typemap)

if q.calculate_residence:
    setupbigtrackingresL(q)
if q.needCEC:
        (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,currenttime,q)                  
fixpbc(tempstate,q)
q.statesL.append(tempstate)

if q.IR: computeIR(q,currenttime)


while True:
    try:
        currenttime=trajtime(trajf,q)  
        tempstate=inputstate(trajf,q,currenttime)  
        if q.needCEC:
            if  q.statesL[-1][0]<currenttime:
                etime=evbtime(evbf)
                (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,etime,q)                                          
                print("%d     %d"%(q.statesL[-1][0],currenttime))
                while etime<currenttime:
                    etime=evbtime(evbf) 
                    (q.CECorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.CECorderL,etime,q)
                    
                                     
        if q.statesL[-1][0]<currenttime:                             
            fixpbc(tempstate,q)
            q.statesL.append(tempstate)    
            if q.IR: computeIR(q,currenttime)
         
                
            
                         

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


        compute(q) #calculate quantities from statesL        
        
       
        kobclean(q.statesL,currenttime)  #clean up states
        
    except StopIteration:           
        break


    

#for i in range(len(q.counthopsL)):
#    print("%d %d\n"%(q.counthopsL[i][0],q.counthopsL[i][1]))

printresults(args[0],q)

msg = 'Enter anything to exit '
uin = raw_input(msg)

