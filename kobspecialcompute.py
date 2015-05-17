import re
import sys
import math
import koblib
from kobcompute import *
import vecmath as vec
from datastruct import *




def construct_resT_CF(q, resL):
    L = []
    for i in range(len(resL)):
        L.append([i * q.smalleststep, 0])
    L[0][1] = resL[0][1]
    for i in range(1, len(resL)):
        for j in range(i, -1, -1):
            L[j][1] += resL[i][1]
    N = L[0][1]
    for item in L:
        item[1] = item[1] * 1. / N
    return L


def calc_meanresT(resL):
    Normalization = 0
    sum = 0
    for item in resL:
        sum += item[0] * item[1]
        Normalization += item[1]
    return sum * 1. / Normalization


def calcResT(q, trackingresL, resL):
    pair = resL[0]
    list1 = q.Typemap[pair[0][0]]
    list2 = q.Typemap[pair[1][0]]
    maxdist = pair[2]
    t_star = pair[3]
    currentt = q.statesL[-1][0]
    for item in pair[0][1:]:
        list1 += q.Typemap[item]
    for item in pair[1][1:]:
        list2 += q.Typemap[item]

    pos1 = 0
    for i in list1:
        pos2 = 0
        for j in list2:
            if j != i:
                xdiff1 = minimagedist(q.statesL[-1][j][2], q.statesL[-1][i][2], q.xboxlength)
                ydiff1 = minimagedist(q.statesL[-1][j][3], q.statesL[-1][i][3], q.yboxlength)
                zdiff1 = minimagedist(q.statesL[-1][j][4], q.statesL[-1][i][4], q.zboxlength)
                # if q.statesL[-1][i][0]==1341 and q.statesL[-1][j][0]==5878:
                #    if dist1<maxdist:
                #        q.debugrestf2.write(" %d  , %d , %f , %d, YES\n "%(q.statesL[-1][i][0],q.statesL[-1][j][0],dist1, currentt,))
                #    else:
                #        q.debugrestf2.write(" %d  , %d , %f , %d, NO\n "%(q.statesL[-1][i][0],q.statesL[-1][j][0],dist1, currentt,))
                dist1 = math.sqrt(xdiff1 ** 2 + ydiff1 ** 2 + zdiff1 ** 2)
                startt = trackingresL[pos1][pos2][0]
                lastt = trackingresL[pos1][pos2][1]
                if dist1 < maxdist:
                    if startt == -1:
                        trackingresL[pos1][pos2][0] = currentt  #initial time in residence
                        trackingresL[pos1][pos2][1] = currentt  #last time in residence
                    else:
                        if currentt - (
                            lastt + q.smalleststep) < t_star:  # note: dist1<maxdist already confirms there is residence in this step.
                            trackingresL[pos1][pos2][1] = currentt  #update the "last" time.
                else:
                    if lastt != -1 and currentt - lastt >= t_star:
                        # note: dist1>=maxdist already confirms there is no residence in this time.
                        #================================================================
                        if startt > q.statesL[0][
                            0] + t_star:  # remember that the first time this script actually reads is the second time slot.
                            # this also prevents the counting the very first residence that could be only partially capatured.
                            timediff = (lastt + q.smalleststep) - startt
                            maxbin = len(resL) - 1
                            whichbin = timediff / q.smalleststep
                            if whichbin > maxbin:
                                while maxbin < whichbin:  # -1 is for the label for the pair at the first slot
                                    if len(resL) == 1:
                                        resL.append([timediff, 0])
                                    else:
                                        resL.append([resL[-1][0] + q.smalleststep, 0])
                                    maxbin += 1
                                resL[-1][1] += 1  #last entry
                            else:
                                resL[whichbin][1] += 1
                                #q.debugrestf.write("%d, %d , %d, %d , %d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],\
                                #   trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep))                            
                        trackingresL[pos1][pos2][0] = -1
                        trackingresL[pos1][pos2][1] = -1
            pos2 += 1
        pos1 += 1


def calcwaterResT(q, trackingresL, resL):
    list = q.Typemap[q.OWID]
    # list=[116,118,175,413,485,756] #for debugging with short.lammpstrj in 80CsFCl
    firstHindex = q.Typemap[q.HWID][0]
    maxdist = q.WresOOdist
    Hbondangle = q.WresHOOangle
    t_star = q.Wreststar
    currentt = q.statesL[-1][0]
    
    numberofHbonds = 0
    resL[0] += 1
    OOsqdist2L=[]
    BlockerOsqdist2L=[]
    
    calcsqdistmatrixsametype(OOsqdist2L,q,[q.OWID])    
    calcsqdistmatrixtwotype(BlockerOsqdist2L,q,[q.OWID],q.blockerL[0])
    for pos1 in range(len(list)):
        i=list[pos1]    
        O1coords = q.statesL[-1][i][q.xl:q.xl + 3]
        #H1O1coords = q.statesL[-1][firstHindex + 2 * pos1][q.xl:q.xl + 3]
        #H2O1coords = q.statesL[-1][firstHindex + 2 * pos1 + 1][q.xl:q.xl + 3]


        for pos2 in range(len(list)):
            j=list[pos2]
            gotHbondangle = False            
            if j > i:
                startt = trackingresL[pos1][pos2][0]
                lastt = trackingresL[pos1][pos2][1]
                O2coords = q.statesL[-1][j][q.xl:q.xl + 3]
                #H1O2coords = q.statesL[-1][firstHindex + 2 * pos2][q.xl:q.xl + 3]
                #H2O2coords = q.statesL[-1][firstHindex + 2 * pos2 + 1][q.xl:q.xl + 3]
                #sqdist1=minimage3Dsqdist(O1coords,O2coords,q.xboxlength)
#                print pos1, pos2
                sqdist1 = OOsqdist2L[pos1][pos2]             
                if sqdist1 < maxdist*maxdist:
                    if q.collectblockerdata:
                        blockersqmindistindex = BlockerOsqdist2L[pos1][-1]
                        blockersqmindist= BlockerOsqdist2L[pos1][blockersqmindistindex]                        
                        q.h.adddata(math.sqrt(blockersqmindist),math.sqrt(sqdist1)) 
                    #                    if vec.angle(H1O1coords,O1coords,O2coords,q.xboxlength)*180/math.pi < Hbondangle:
                    #                         gotHbondangle=True
                    #                    elif vec.angle(H2O1coords,O1coords,O2coords,q.xboxlength)*180/math.pi< Hbondangle:
                    #                         gotHbondangle=True
                    #                    elif vec.angle(H1O2coords,O2coords,O1coords,q.xboxlength)*180/math.pi< Hbondangle:
                    #                         gotHbondangle=True
                    #                    elif vec.angle(H2O2coords,O2coords,O1coords,q.xboxlength)*180/math.pi< Hbondangle:
                    #                         gotHbondangle=True
                         
                    #if dist1 < maxdist and gotHbondangle:
                    numberofHbonds += 2
                    if startt == -1:
                        #                             unsorteddataL=[data1,data2,data3,data4]
                        #                             dataL=sorted(unsorteddataL)
                        #                             minangle=dataL[0][0]
                        #                             Hindex=dataL[0][1]
                        trackingresL[pos1][pos2][0] = currentt
                        trackingresL[pos1][pos2][1] = currentt
                        trackingresL[pos1][pos2][2] = -999
                    else:
                        if currentt - (
                            lastt + q.smalleststep) < t_star:  # note: dist1<maxdist already confirms there is residence in this step.
                            trackingresL[pos1][pos2][1] = currentt
                        #trackingresL[pos1][pos2][3].append(dist1)
                else:
                    if lastt != -1 and currentt - lastt >= t_star:
                        # note: dist1>=maxdist already confirms there is no
                        # residence in this time.
                        #================================================================
                        if startt > q.statesL[0][
                            0] + t_star:  # remember that the first time this script actually reads is the second time
                            # slot.
                            if sqdist1 < maxdist*maxdist:
                                q.lostanglecount += 1  #only angle broke

                            timediff = (lastt + q.smalleststep) - startt

                            maxbin = len(resL) - 1
                            whichbin = timediff / q.Hbondbinsize
                            if whichbin > maxbin:
                                while maxbin < whichbin:  # -1 is for the label for the pair at the first slot
                                    if len(resL) == 1:
                                        resL.append([timediff, 0])
                                    else:
                                        resL.append([resL[-1][0] + q.smalleststep, 0])
                                    maxbin += 1
                                resL[-1][1] += 1  #last entry
                            else:
                                resL[whichbin][1] += 1
                            if q.collectblockerdata:

                                targetdataL = q.blockerdataL[0]
                                L = trackingresL[pos1][pos2][3]
                                for item in L:
                                    whichbin = int(item / 0.1) + 1
                                    koblib.histo(targetdataL, 0.1, whichbin, 1)
                                #                                    targetdataL=q.blockerdataL[1]
                                #                                    lastdist=trackingresL[pos1][pos2][3][-1]
                                #                                    whichbin=int(lastdist/0.1)+1
                                #                                    koblib.histo(targetdataL,0.1,whichbin,1)

                        trackingresL[pos1][pos2][0] = -1
                        trackingresL[pos1][pos2][1] = -1
                        trackingresL[pos1][pos2][2] = None
                        trackingresL[pos1][pos2][3] = []
    q.avgnoHbonds += numberofHbonds * 1. / len(list)


def calcFPTRDF(q, RDFL):
    type1 = RDFL[0][0][0]
    type2 = RDFL[0][1][0]
    list1 = []
    list2 = None
    currtime = q.statesL[-1][0]
    starttime = q.statesL[0][0]
    whichtimeseriesindex = -999
    whichtimebin = (currtime - starttime) / q.smalleststep - 1
    count = 0

    for item in q.FPTtypes:
        type = item[0]
        if type1 == type:
            whichtimeseriesindex = count
        count += 1
    if whichtimeseriesindex == -999:
        sys.exit("can't find whichtimeseriesindex")

    list1 = q.timeseries[whichtimeseriesindex][whichtimebin][1]
    RDFL[0][2] += len(q.timeseries[whichtimeseriesindex][whichtimebin][1])

    list2 = q.Typemap[type2]

    for i in list1:
        # print i
        atominfo_i = q.statesL[-1][i]
        for j in list2:
            if j != i:
                xdiff = minimagedist(q.statesL[-1][j][2], atominfo_i[2], q.xboxlength)
                ydiff = minimagedist(q.statesL[-1][j][3], atominfo_i[3], q.yboxlength)
                zdiff = minimagedist(q.statesL[-1][j][4], atominfo_i[4], q.zboxlength)
                dist = math.sqrt(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)

                whichbin = int(dist / q.binsize) + 1  # +1 to account for the label slot
                koblib.histo(RDFL, q.binsize, whichbin, 1)


def findFPT(q):
    count = 0
    currtime = q.statesL[-1][0]
    firstt = q.statesL[0][0]

    for item in q.FPTtypes:
        MFPT = item[2]
        if MFPT != -999:
            q.timeseries[count].append([currtime, []])
        type = item[0]
        threshold = item[1]
        if len(q.bigtrackingFPTL[count]) != len(q.Typemap[type]) + 1:  # +1 for label
            for item in q.Typemap[type]:
                q.bigtrackingFPTL[count].append([currtime, q.statesL[-1][item][q.xl:q.xl + 3]])
        else:
            index = 1
            for item in q.Typemap[type]:
                currpos = q.statesL[-1][item][q.xl:q.xl + 3]
                oldtime = q.bigtrackingFPTL[count][index][0]
                oldpos = q.bigtrackingFPTL[count][index][1]
                dist = minimage3Ddist(currpos, oldpos, q.boxlengthL)
                if dist >= threshold:
                    timediff = currtime - oldtime
                    whichbin = timediff / q.FPTbinsize + 1
                    koblib.histo(q.bigFPTL[count + 1], q.FPTbinsize, whichbin, 1)
                    q.bigtrackingFPTL[count][index][0] = currtime
                    q.bigtrackingFPTL[count][index][1] = currpos
                    if MFPT != -999:
                        if timediff >= MFPT:
                            whichtimeslot = (oldtime - firstt) / q.smalleststep - 1
                            howmanyslots = timediff / q.smalleststep + 1
                            for i in range(howmanyslots):
                                q.timeseries[count][whichtimeslot + i][1].append(item)
                index += 1
        count += 1


def specialcompute(q):
    if q.findFPT:
        currtime = q.statesL[-1][0]
        q.bigFPTRDFL[0] += 1
        count = 1  # the 0th place is counter
        for pair in q.FPTpairL:
            calcFPTRDF(q, q.bigFPTRDFL[count])
            count += 1
    if q.calculate_residence:
        count = 0
        for pair in q.respairL:
            calcResT(q, q.bigtrackingresL[count], q.bigresL[count])
            count += 1
    if q.calculate_waterresidence:
        calcwaterResT(q, q.watertrackingresL, q.waterresL)

    if q.findFPT:
        findFPT(q)

    if q.calculate1DSqdisp:
        typecount=1
        for type in q.Sqdisp:
            fout=open('1DSq_'+str(type)+'_'+suffix+'.csv','w')
            maxlen=0
            itemcount=1
            index=1
            for item in q.big1DSqdisp[typecount][1:]:
                if len(item[1])-1>maxlen:
                    maxlen=len(item[1])
                    index=itemcount
                itemcount+=1
            fout.write("r ,")
            for item in q.big1DSqdisp[typecount][1:]:
                fout.write(" t=%d, "%item[0])
            fout.write("\n")
            for i in range(1,maxlen):
                fout.write(" %f ,"%q.big1DSqdisp[typecount][index][1][i][0])
                for item in q.big1DSqdisp[typecount][1:]:
                    if i < len(item[1]):                        
                        actualcount=item[1][i][1]*1.                        
                        fout.write(" %f ,"%(actualcount))
                    else:
                        fout.write(" -999 ,")
                fout.write("\n")            
            typecount+=1
            fout.close()    
