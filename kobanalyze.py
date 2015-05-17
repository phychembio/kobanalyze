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


OWID=1
q.OWID=OWID
HWID=2
q.HWID=HWID
OHID=3
q.OHID=OHID
HHID=4
q.HHID=HHID

prefix="hydroxide/"

q.debugf=open("reactimes.txt","w")  


q.trajfile=prefix+'out'+args[0]
q.evbfile=prefix+'out'+args[0]
q.molfile=prefix+'start.lammpstrj'
driftfile='drift.data'
q.findFPT=False
q.calculate_waterresidence=False
q.calculate_residence=False
q.calculate_interface=False
q.calc_reactivewater=True

q.needCEC=True
q.MSDtypes=[1]
q.pairL=[]
#q.pairL=[([NID],[ClID]),([NID],[FID]),([NID],[OWID]),([FID],[ClID]),([FID],[OWID]),([ClID],[OWID])]                                                      


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



if q.needCEC:
        (q.CECorderL)=insertCEC(tempstate,evbf,q.OorderL,q.CECorderL,currenttime,q)


                   
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
                (q.CECorderL)=insertCEC(tempstate,evbf,q.OorderL,q.CECorderL,etime,q)                          
                #print("%d     %d"%(q.statesL[-1][0],currenttime))
                while etime<currenttime:
                    etime=evbtime(evbf) 
                    (q.CECorderL)=insertCEC(tempstate,evbf,q.OorderL,q.CECorderL,etime,q)
                                     
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


