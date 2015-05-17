#!/usr/bin/python

#modified on 5/23/2013
#added modules koblib and kobcompute

import re
import math
import sys
import numpy
from koblib import *
from kobcompute import *
from kobinfo import*


from operator import itemgetter

OWID=7
HWID=8
FID=9
NID=6


trajfile='3p/3HClshort'
evbfile='3p/3HClshort'
#trajfile="shortewald"
#evbfile="shortewald"
trajext='.lammpstrj'
evbext='.evb'
molfile='start.lammpstrj'
restartfile='restart.kob'
smalleststep=findprintfreq(trajfile+trajext)
binsize=0.1
noofCEC=None
skipstep=0
needmolID=False #needed for reactive dynamics
computeC1Sdist=False
restart=False
MSDdecomposition=False
calculateclosest=False
c_i_info=False;

q=KobInfo(trajfile,evbfile,trajext,evbext,molfile,restartfile,smalleststep,binsize,noofCEC,skipstep,needmolID,computeC1Sdist,restart,MSDdecomposition,calculateclosest,c_i_info)
#q=KobInfo('shortewald','shortewald','.lammpstrj','.evb','start.lammpstrj','restart.kob',1000,0.1,40,0,False,False,False,False,False,True)


q.MSDtypes=[1,-1]
q.pairL=[([1],[-1])]
q.closestL=[([FID],[NID]),([NID],[999])]
q.respairL=[]
q.notpolymerL=[OWID,HWID,NID,98] #for calculating "closest", these are not included

initialize(q)

stepskipped=0 #time steps already skipped, no need to change

args=sys.argv[1:]
if len(args)==0:
    args.append("")


trajf=open(q.trajfile+args[0]+q.trajext,'rU')         
if q.needCEC:
    evbf=open(q.evbfile+args[0]+q.evbext,'rU') 
else:
    evbf=None
    
          
    
if q.needCEC:        
    currenttime=eventime(evbf,trajf,q)
    etime=currenttime      
else:
    currenttime=trajtime(trajf,q)  
####End initialize#####    

tempstate=inputstate(trajf,currenttime)
if q.needCEC:
        (q.orderL,q.lastorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.orderL,currenttime,q)
                   
fixpbc(tempstate,q.statesL)
q.statesL.append(tempstate)


setupTypemap(q.statesL,q.Typemap)


while True:
    try:
        currenttime=trajtime(trajf,q)  
        tempstate=inputstate(trajf,currenttime)  
        if q.needCEC:
            if  q.statesL[-1][0]<currenttime:
                etime=evbtime(evbf)
                (q.orderL,q.lastorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.orderL,etime,q)                          
                print("%d     %d"%(q.statesL[-1][0],currenttime))
                while etime<currenttime:
                    etime=evbtime(evbf) 
                    (q.orderL,q.lastorderL)=insertCEC(tempstate,evbf,q.waterorderL,q.orderL,etime,q)
                                     
        if q.statesL[-1][0]<currenttime:                             
            fixpbc(tempstate,q.statesL)
            q.statesL.append(tempstate)             

####Check/Debug########
        if currenttime%10000==0:
            print currenttime
            sys.stdout.flush()            
        if q.MSDdecomposition==False and len(q.statesL)>0 and len(q.statesL[-1])!=len(q.statesL[0])\
            or q.MSDdecomposition==True and len(q.statesL)>2 and len(q.statesL[-1])+2*noofCEC!=len(q.statesL[0]):
            print "Current Time: " +str(currenttime)
            print "latest state size: " +str(len(q.statesL[-1])) 
            print "initial state size: " +str(len(q.statesL[0]))
            sys.exit("state sizes mismatch..quitting")     
           
####End Check/Debug########


        compute(q) #calculate quantities from statesL

        kobclean(q.statesL,currenttime)  #clean up states
        
    except StopIteration:           
        break


    



printresults(args[0],q)

msg = 'Enter anything to exit '
uin = raw_input(msg)

