#!/usr/bin/env python

#modified on 5/23/2013
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

q=KobInfo()

OWID=7
HWID=2
OHID=8
FID=9
NID=6
ClID=10

prefix="PVBTMA-FCl/"
q.trajfile=prefix+'short'+args[0]
q.evbfile=prefix+'short'+args[0]
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
q.c_i_info=False
q.calculate_residence=False
q.SID=10
q.C1ID=2
q.OID=5

#q.debugrestf=open("debug_shortFCl_N-Cl3.csv","w")
#q.debugrestf2=open("debug_shortFCl_1023-5891.csv","w")
#q.debugrestf.write("ps, N_ID, Cl_ID, start time, end time\n")
#q.debugrestf2.write("N_ID, Cl_ID, dist, timestep, smaller than 6.5?\n")
q.MSDtypes=[2]
q.pairL=[([2],[2])]
#q.pairL=[([NID],[ClID]),([NID],[FID]),([NID],[OWID]),([FID],[ClID]),([FID],[OWID]),([ClID],[OWID])]                                                      
if q.calculate_residence:
    q.respairL=[([NID],[ClID],6.5,2000),([NID],[FID],6.5,2000),([NID],[FID],8.3,2000)] #([A],[B],Max R_A-B, t*)                                              
#q.respairL=[([NID],[FID],7.25,2000)] #([A],[B],Max R_A-B, t*)  
if q.calculateclosest:
    q.closestL=[([FID],[NID], True),([ClID],[NID], True),([NID],[NID], False)]
    q.conditionL=[([FID],[OWID],3.25),([ClID],[OWID],3.9)]
    q.notpolymerL=[OWID,HWID,NID,98] #for calculating "closest", these are not included


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
#        q.Typemap[:]=[[]]
#        setupTypemap(tempstate,q.Typemap) #not a repeat to the first call. Some info from the first Typemap is needed before the CECs are inserted.
                   
fixpbc(tempstate,q)
q.statesL.append(tempstate)


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

####Check/Debug########
        if currenttime%10000==0:
            print currenttime
            sys.stdout.flush()            
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

