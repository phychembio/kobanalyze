import re
import sys
import math
import koblib
from kobcompute import *
from kobspecialcompute import *
import vecmath as vec
from datastruct import *


def setupbigtrackingresL(q):
    count = 0
    for resL in q.bigresL:
        N_A = len(q.Typemap[resL[0][0][0]])
        N_B = len(q.Typemap[resL[0][1][0]])
        for i in range(N_A):
            q.bigtrackingresL[count].append([])
            for j in range(N_B):
                q.bigtrackingresL[count][-1].append([-1, -1])
        count += 1


def setupwatertrackingresL(q):
    N_A = len(q.Typemap[q.OWID])
    q.lostanglecount = 0
    for i in range(N_A):
        q.watertrackingresL.append([])
        for j in range(N_A):
            q.watertrackingresL[-1].append([-1, -1, -1, []])
    if q.collectblockerdata:
        q.blockerindexL = []
        q.blockerdataL = [[q.blockerL[0]], [q.blockerL[0]]]
        for item in q.blockerL[0]:
            q.blockerindexL += q.Typemap[item]

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

def specialinit(q):
    if q.calculate_residence:
        q.bigresL = []
        q.bigtrackingresL = []
        for pair in q.respairL:
            q.bigresL.append([pair])
            if -1 in pair[0] or -1 in pair[1]:
                q.needCEC = True
            q.bigtrackingresL.append([])
        setupbigtrackingresL(q)

    if q.calculate_waterresidence:
        q.waterresL = [0]
        q.watertrackingresL = []
        q.avgnoHbonds = 0
        q.h = hist2d()
        setupwatertrackingresL(q)

    q.bigFPTRDFL = []
    q.bigFPTRDFL.append(0)
    q.stopcalcuteFPT = True
    for pair in q.FPTpairL:
        pair = list(pair)
        pair.append(0)
        q.bigFPTRDFL.append([pair])

    if q.calculateFPTproperties:
        count = 0
        for item in q.FPTtypes:
            type = item[0]
            MFPT = item[2]
            inf = open('Slow-%d_type=' % MFPT + str(type) + '.txt', 'r')
            for line in inf:
                info = line.split()
                time = int(info[0])
                atomL = [int(x) for x in info[1:]]
                q.timeseries[count].append([time, atomL])
            count += 1
        q.findFPT = False
    if q.findFPT:
        q.timeseries = []
        q.bigFPTL = [0]
        for item in q.FPTtypes:
            q.timeseries.append([])
            type = item[0]
            q.bigFPTL.append([type])
            q.bigtrackingFPTL.append([type])
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


def specialprintresults(suffix, q):
    if q.findFPT:  # writing out long stays
        count = 0
        for item in q.FPTtypes:
            type = item[0]
            MFPT = item[2]
            FPT = q.bigFPTL[count + 1]
            FPTf = open('FPT_type=' + str(type) + suffix + '.txt', 'w')
            sum = 0
            for bin in FPT[1:]:
                sum += bin[1]
                FPTf.write("%d  %d \n" % (bin[0] + q.FPTbinsize / 2, bin[1]))
            FPTf.close()
            if MFPT != -999:
                FPTf = open('Slow-%d_type=' % MFPT + str(type) + '.txt', 'w')
                L = q.timeseries[count]
                for item in L:
                    FPTf.write("%d %s\n" % (item[0], " ".join(map(str, item[1]))))
                FPTf.close()
            count += 1

    if q.findFPT:  # RDF
        count = 1  #the 0th place is counter
        volume = q.xboxlength * q.yboxlength * q.zboxlength
        for pair in q.FPTpairL:
            nooftype2 = countatoms(q.Typemap, pair[1])
            RDF = q.bigFPTRDFL[count]
            rdff = open('FPTRDF_type=' + str(pair[0]) + '_' + str(pair[1]) + suffix + '.txt', 'w')
            sum = 0
            type1count = RDF[0][2]
            for bin in RDF[1:]:
                sum += bin[1]
                intensity = bin[1] / (4 * math.pi * bin[0] ** 2 * q.binsize * nooftype2 / volume * type1count)
                rdff.write("%f   %f\n" % (bin[0], intensity))
            rdff.close()
            count += 1

    if q.calculate_residence:
        for i in range(len(q.respairL)):
            resf = open('ResTCF_type=' + str(q.respairL[i][0]) + '-' + str(q.respairL[i][1]) + '-' + str(
                q.respairL[i][2]) + '_' + suffix + '.txt', 'w')
            # resCF=kobcompute.construct_resT_CF(q,q.bigresL[i][1:])
            resCF = q.bigresL[i][1:]
            meanresT = calc_meanresT(q.bigresL[i][1:])
            resf.write("Mean_Residence_Time= %f\n" % meanresT)
            for item in resCF:
                resf.write("%d   %f \n" % ( item[0], item[1]))
            resf.close()

    if q.calculate_waterresidence:
        resf = open('HbondResTCF_' + suffix + '.txt', 'w')
        # resCF=kobcompute.construct_resT_CF(q,q.bigresL[i][1:])
        resCF = q.waterresL[1:]
        totalcount = 0
        for item in resCF:
            resf.write("%d   %f \n" % ( item[0], item[1]))
            totalcount += item[1]
        resf.write("Total count: %d\n" % totalcount)
        resf.write("Avg no of Hbonds: %f\n" % (q.avgnoHbonds / q.waterresL[0]))

        meanresT = calc_meanresT(q.waterresL[1:])
        resf.write("Mean_Residence_Time= %f\n" % meanresT)
        resf.close()

        if q.collectblockerdata:
            q.h.output("ClosestOOClosestBlockerO"+suffix+".csv")
            # blockerresf = open('HbondBlockerdist-OOdist_' + suffix + '.txt', 'w')
            # totalcount = 0
            # for item in q.blockerdataL[0][1:]:
            #     blockerresf.write("%f   %f \n" % ( item[0], item[1]))
            #     totalcount += item[1]
            # blockerresf.write("Total count: %d\n" % totalcount)
            # blockerresf.write("Average Hbond lifetimeset: %f\n" % q.Hbondlife)
            # blockerresf.close()
            # # blockerresf=open('HbondBlockerdist-last_'+suffix+'.txt','w')
            # # totalcount=0
            # # for item in q.blockerdataL[1][1:]:
            # # blockerresf.write("%f   %f \n" %( item[0] , item[1]))
            # #     totalcount+=item[1]
            # # blockerresf.write("Total count: %d\n"%totalcount)
            # # blockerresf.write("Average Hbond lifetimeset: %f\n"%q.Hbondlife)
            # blockerresf.close()
