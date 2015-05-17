import re
import sys
import math
import koblib
from kobcompute import *
from kobspecialcompute import *
import vecmath as vec



def setupbigtrackingresL(q):
    count=0
    for resL in q.bigresL:
        N_A=len(q.Typemap[resL[0][0][0]])
        N_B=len(q.Typemap[resL[0][1][0]])
        for i in range(N_A):
            q.bigtrackingresL[count].append([])
            for j in range(N_B):
                q.bigtrackingresL[count][-1].append([-1,-1])
        count+=1

def setupwatertrackingresL(q):    
    N_A=len(q.Typemap[q.OWID])
    q.lostanglecount=0    
    for i in range(N_A):
        q.watertrackingresL.append([])
        for j in range(N_A):
            q.watertrackingresL[-1].append([-1,-1,-1]) 
    if q.collectblockerdata:
        q.blockerindexL=[]
        q.blockerdataL=[q.blockerL[0]]
        for item in q.blockerL[0]:
            q.blockerindexL+=q.Typemap[item]
        


def specialinit(q):
    if q.calculate_residence:
        q.bigresL=[]
        q.bigtrackingresL=[]
        for pair in q.respairL:
            q.bigresL.append([pair])         
            if -1 in pair[0] or -1 in pair[1]:
                q.needCEC=True
            q.bigtrackingresL.append([])
        setupbigtrackingresL(q)

    if q.calculate_waterresidence:
        q.waterresL=[[q.OWID,q.OWID]]
        q.watertrackingresL=[]
        setupwatertrackingresL(q)
        
    q.bigFPTRDFL=[]
    q.bigFPTRDFL.append(0)
    q.stopcalcuteFPT=True
    for pair in q.FPTpairL:
        pair=list(pair)
        pair.append(0)        
        q.bigFPTRDFL.append([pair])
        
    if q.calculateFPTproperties:
        count=0
        for item in q.FPTtypes:
            type=item[0]
            MFPT=item[2]            
            inf=open('Slow-%d_type='%MFPT+str(type)+'.txt','r')            
            for line in inf:
                info=line.split()
                time=int(info[0])
                atomL=[int(x) for x in info[1:]]
                q.timeseries[count].append([time,atomL])
            count+=1
        q.findFPT=False 
           

def specialprintresults(suffix,q):

    if q.findFPT:
        count=1  #the 0th place is counter    
        volume=q.xboxlength*q.yboxlength*q.zboxlength
        for pair in q.FPTpairL:                
            nooftype2=countatoms(q.Typemap,pair[1])
            RDF=q.bigFPTRDFL[count]
            rdff=open('FPTRDF_type='+str(pair[0])+'_'+str(pair[1])+suffix+'.txt','w')
            sum=0
            type1count=RDF[0][2]
            for bin in RDF[1:]:
                sum+=bin[1]
                intensity=bin[1]/(4*math.pi*bin[0]**2*q.binsize*nooftype2/volume*type1count)                          
                rdff.write("%f   %f\n" %(bin[0] ,intensity))
            rdff.close() 
            count+=1

    if q.calculate_residence:
        for i in range(len(q.respairL)):
            resf=open('ResTCF_type='+str(q.respairL[i][0])+'-'+str(q.respairL[i][1])+'-'+str(q.respairL[i][2])+'_'+suffix+'.txt','w')
           # resCF=kobcompute.construct_resT_CF(q,q.bigresL[i][1:])
            resCF=q.bigresL[i][1:]
            meanresT=kobcompute.calc_meanresT(q.bigresL[i][1:])
            resf.write("Mean_Residence_Time= %f\n"%meanresT)
            for item in resCF:            
                resf.write("%d   %f \n" %( item[0] , item[1]))
            resf.close()
    if q.calculate_waterresidence:        
        resf=open('WaterResTCF_'+suffix+'.txt','w')
        # resCF=kobcompute.construct_resT_CF(q,q.bigresL[i][1:])
        resCF=q.waterresL[1:]        
        totalcount=0
        for item in resCF:            
            resf.write("%d   %f \n" %( item[0] , item[1]))
            totalcount+=item[1]
        resf.write("Total count: %d\n"%totalcount)
        resf.write("Lost angle count: %d\n"%q.lostanglecount)
        meanresT=calc_meanresT(q.waterresL[1:])
        resf.write("Mean_Residence_Time= %f\n"%meanresT)
        resf.close()
        if q.collectblockerdata:
            blockerresf=open('WaterBlockerdist_'+suffix+'.txt','w')
            totalcount=0
            for item in q.blockerdataL[1:]:            
                blockerresf.write("%f   %f \n" %( item[0] , item[1]))
                totalcount+=item[1]
            blockerresf.write("Total count: %d\n"%totalcount)