import re
import sys
import math
import koblib
from kobcompute import *
import vecmath as vec


def findzcenter1(q): #gibbs dividing surface across all xy-plane
     noofatom=q.statesL[-1][1]
     binsize=0.5
     bulkdensity=0.63
     binvol=q.xboxlength*q.yboxlength*binsize
     

     zhisto=["local"]
     for i in range(noofatom):
        z=q.statesL[-1][3+i][q.xl+2]
        if z>0:
            typeIndex=q.statesL[-1][3+i][1]-1
            whichbin=int(z/binsize)+1
            koblib.histo(zhisto,binsize,whichbin,q.massL[typeIndex]/binvol)    
     firstboundary=-1
     secondboundary=-1
     count=1
     for item in zhisto[1:]:
         density=0
         if item[2]>0:
            density=item[1]
            if firstboundary==-1:
                if density>=bulkdensity/2:
                    deltadensity=bulkdensity-item[1]
                    slope=(zhisto[count][1]-zhisto[count-1][1])/(zhisto[count][0]-zhisto[count-1][0])
                    firstboundary=item[0]-slope*deltadensity
            elif secondboundary==-1:
                if density<=bulkdensity/2 and item[0]-firstboundary>5:
                    deltadensity=bulkdensity-item[1]
                    slope=(zhisto[count][1]-zhisto[count-1][1])/(zhisto[count][0]-zhisto[count-1][0])
                    secondboundary=item[0]+slope*deltadensity
         count+=1
     


     return (secondboundary+firstboundary)/2

def findzcenter2(q): #using z center of mass
     noofatom=q.statesL[-1][1]
     sum=0
     masssum=0
     for i in range(noofatom):
         type=q.statesL[-1][i+3][1]
         zpos=q.statesL[-1][i+3][4]         
         mass=q.massL[type-1]
         if q.initialzcenterofmass!=-1:
            if abs(zpos-q.initialzcenterofmass)<50: #not letting gas molecules affect this
                 sum+=q.massL[type-1]*zpos
                 masssum+=mass
         else:
            sum+=q.massL[type-1]*zpos
            masssum+=mass
     if(q.initialzcenterofmass)==-1:
         q.initialzcenterofmass=sum/masssum     
     return sum/masssum

def findzcenter3(q): #using z center of mass for each \Delta x \Delta y
     noofatom=q.statesL[-1][1]
     
     xwidth=q.xhi-q.xlo
     ywidth=q.yhi-q.ylo
     noofxybins=10
     xbinsize=xwidth/noofxybins
     ybinsize=ywidth/noofxybins

     sum=[]
     for i in range(noofxybins):
         L=[]
         for j in range(noofxybins):
             L.append(0)
         sum.append(L)

     masssum=[]
     for i in range(noofxybins):
         L=[]
         for j in range(noofxybins):
             L.append(0)
         masssum.append(L)
     zposL=[]
     for i in range(noofxybins):
         L=[]
         for j in range(noofxybins):
             L.append(0)
         zposL.append(L)

     for i in range(noofatom):
         type=q.statesL[-1][i+3][1]
         xpos=q.statesL[-1][i+3][2]         
         if xpos<=q.xlo:
             while xpos<=q.xlo:
                 xpos+=xwidth
         elif xpos>=q.xhi:
             while xpos>=q.xhi:
                 xpos-=xwidth
         ypos=q.statesL[-1][i+3][3]         
         if ypos<=q.ylo:
             while ypos<=q.ylo:
                 ypos+=ywidth
         elif ypos>=q.yhi:
             while ypos>=q.yhi:
                 ypos-=ywidth
         whichxbin=int(xpos/xbinsize)
         whichybin=int(ypos/ybinsize)
         zpos=q.statesL[-1][i+3][4]         
         mass=q.massL[type-1]                  
     
         if(q.initialzcenterofmass)!=-1:
             if abs(zpos-q.initialzcenterofmass)<50:
                 sum[whichxbin][whichybin]+=mass*zpos
                 masssum[whichxbin][whichybin]+=mass
         else:
             sum[whichxbin][whichybin]+=mass*zpos
             masssum[whichxbin][whichybin]+=mass

     if(q.initialzcenterofmass)==-1:
         tempsum=0
         tempmasssum=0
         for i in range(noofxybins):
            for j in range(noofxybins):
                tempsum+=sum[i][j]
                tempmasssum+=masssum[i][j]
         q.initialzcenterofmass=tempsum/tempmasssum
     for i in range(noofxybins):
        for j in range(noofxybins):
            zposL[i][j]=sum[i][j]/masssum[i][j]

     return zposL


     return (secondboundary+firstboundary)/2

def interfacelocalcompute(q):
     
     noofatom=q.statesL[-1][1]
     
     variance=0
     xwidth=q.xhi-q.xlo
     ywidth=q.yhi-q.ylo
     noofxybins=10
     xbinsize=xwidth/noofxybins
     ybinsize=ywidth/noofxybins
     
     zcenterL=findzcenter3(q)

     q.bigionL[0]+=1
     ionbinsize=0.1

     for type in q.iontype:
         for index in q.Typemap[type]:
             if type!=4 or (type==4 and (index-q.Typemap[4])%6==0):                              
                 for L in q.bigionL[1:]:
                     if type==L[0]:
                         xpos=q.statesL[-1][index][2]         
                         if xpos<=q.xlo:
                             while xpos<=q.xlo:
                                 xpos+=xwidth
                         elif xpos>=q.xhi:
                             while xpos>=q.xhi:
                                 xpos-=xwidth
                         ypos=q.statesL[-1][index][3]         
                         if ypos<=q.ylo:
                             while ypos<=q.ylo:
                                 ypos+=ywidth
                         elif ypos>=q.yhi:
                             while ypos>=q.yhi:
                                 ypos-=ywidth
                         whichxbin=int(xpos/xbinsize)
                         whichybin=int(ypos/ybinsize)
                         zcenter=zcenterL[whichxbin][whichybin]
                         zpos=q.statesL[-1][index][q.xl+2]+q.statesL[-1][index][q.xl+5]*q.zboxlength
                         #zcenter=zcenterL[whichxbin][whichybin]
                         if zpos-zcenter>0:
                             whichbin=int((zpos-zcenter)/ionbinsize)+1
                         else:
                             whichbin=int((zcenter-zpos)/ionbinsize)+1
                         koblib.histo(L,ionbinsize,whichbin,1)


def interfacecompute(q):
     noofatom=q.statesL[-1][1]     
     variance=0               
     zcenter=findzcenter2(q)

     #zcenter:=findzcenter3(q)

     q.bigionL[0]+=1
     ionbinsize=0.5

     for type in q.iontype:
         for index in q.Typemap[type]:
             if type!=4 or (type==4 and (index-q.Typemap[4][0])%6==0):    
                for L in q.bigionL[1:]:
                    if type==L[0]:
                        zpos=q.statesL[-1][index][q.xl+2]+q.statesL[-1][index][q.xl+5]*q.zboxlength
                        if zpos-zcenter>0:
                            whichbin=int((zpos-zcenter)/ionbinsize)+1
                        else:
                            whichbin=int((zcenter-zpos)/ionbinsize)+1
                        koblib.histo(L,ionbinsize,whichbin,1)

def construct_resT_CF(q,resL):
    L = []
    for i in range(len(resL)):
        L.append([i * q.smalleststep,0])
    L[0][1] = resL[0][1]
    for i in range(1,len(resL)):
        for j in range(i,-1,-1):
            L[j][1]+=resL[i][1]
    N = L[0][1]
    for item in L:
        item[1] = item[1] * 1. / N
    return L

def calc_meanresT(resL):
    Normalization = 0
    sum = 0
    for item in resL:
        sum+=item[0] * item[1]
        Normalization+=item[1]
    return sum * 1. / Normalization

def calcResT(q,trackingresL,resL):
       pair = resL[0]       
       list1 = q.Typemap[pair[0][0]]
       list2 = q.Typemap[pair[1][0]]
       maxdist = pair[2]
       t_star = pair[3]
       currentt = q.statesL[-1][0]
       for item in pair[0][1:]:
          list1+=q.Typemap[item]
       for item in pair[1][1:]:
           list2+=q.Typemap[item]

       pos1 = 0                           
       for i in list1:            
            pos2 = 0
            for j in list2:              
                if j != i:                
                    xdiff1 = minimagedist(q.statesL[-1][j][2],q.statesL[-1][i][2],q.xboxlength)
                    ydiff1 = minimagedist(q.statesL[-1][j][3],q.statesL[-1][i][3],q.yboxlength)
                    zdiff1 = minimagedist(q.statesL[-1][j][4],q.statesL[-1][i][4],q.zboxlength)
                    dist1 = math.sqrt(xdiff1 ** 2 + ydiff1 ** 2 + zdiff1 ** 2)                                        
       
                startt = trackingresL[pos1][pos2][0]   
                lastt = trackingresL[pos1][pos2][1]        
                if dist1 < maxdist:                           
                    if startt == -1:                         
                         trackingresL[pos1][pos2][0] = currentt
                         trackingresL[pos1][pos2][1] = currentt
                    else:                         
                         if currentt - (lastt + q.smalleststep) < t_star: # note: dist1<maxdist already confirms there is residence in this step.
                                                                     # see the
                                                                                                                                              # difference
                                                                                                                                              # about 20
                                                                                                                                              # lines
                                                                                                                                              # later
                            trackingresL[pos1][pos2][1] = currentt
                         else:
#================================================================
                            if startt > q.statesL[0][0] + t_star: # remember that the first time this script actually reads is the second time
                                                                  # slot.
                                timediff = (lastt + q.smalleststep) - startt                 
                                maxbin = len(resL) - 1 
                                whichbin = timediff / q.smalleststep
                                if whichbin > maxbin:                      
                                    while maxbin < whichbin: # -1 is for the label for the pair at the first slot
                                        if len(resL) == 1:
                                            resL.append([timediff,0])
                                        else:
                                            resL.append([resL[-1][0] + q.smalleststep,0])  
                                        maxbin+=1
                                    resL[-1][1]+=1 #last entry
                                else:                       
                                    resL[whichbin][1]+=1   
                                #q.debugrestf.write("%d, %d , %d, %d ,
                                #%d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep)) #================================================================
                            trackingresL[pos1][pos2][0] = currentt
                            trackingresL[pos1][pos2][1] = currentt
                else:                      
                    if lastt != -1 and currentt - lastt >= t_star:        
                        # note: dist1>=maxdist already confirms there is no
                        # residence in this time.
#================================================================
                        if startt > q.statesL[0][0] + t_star: # remember that the first time this script actually reads is the second time
                                                              # slot.
                            timediff = (lastt + q.smalleststep) - startt                                
                            maxbin = len(resL) - 1 
                            whichbin = timediff / q.smalleststep
                            if whichbin > maxbin:                      
                                while maxbin < whichbin: # -1 is for the label for the pair at the first slot
                                    if len(resL) == 1:
                                        resL.append([timediff,0])
                                    else:
                                        resL.append([resL[-1][0] + q.smalleststep,0])  
                                    maxbin+=1
                                resL[-1][1]+=1 #last entry
                            else:                       
                                resL[whichbin][1]+=1   
                            #q.debugrestf.write("%d, %d , %d, %d ,
                            #%d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep)) #================================================================
                        trackingresL[pos1][pos2][0] = -1
                        trackingresL[pos1][pos2][1] = -1
                pos2+=1
            pos1+=1




def calcwaterResT(q,trackingresL,resL):       
       list = q.Typemap[q.OWID]   
       firstHindex=q.Typemap[q.HWID][0]
       maxdist = q.WresOOdist
       Hbondangle=q.WresHOOangle
       t_star =  q.Wreststar
       currentt = q.statesL[-1][0]              
       
       pos1 = 0                           
       for i in list:                  
            pos2 = 0
            O1coords=q.statesL[-1][i][q.xl:q.xl+3]
            H1O1coords=q.statesL[-1][firstHindex+2*pos1][q.xl:q.xl+3]
            H2O1coords=q.statesL[-1][firstHindex+2*pos1+1][q.xl:q.xl+3]
            for j in list:       
                gotHbondangle=False
                dist1=999                                                                                                   
                if j > i:                                    
                    O2coords=q.statesL[-1][j][q.xl:q.xl+3]
                    H1O2coords=q.statesL[-1][firstHindex+2*pos2][q.xl:q.xl+3]
                    H2O2coords=q.statesL[-1][firstHindex+2*pos2+1][q.xl:q.xl+3]                                        
                    dist1 = minimage3Ddist(O1coords,O2coords,q.xboxlength)
                    if dist1 < maxdist:
                        if vec.angle(H1O1coords,O1coords,O2coords,q.xboxlength)*180/math.pi<Hbondangle  or\
                           vec.angle(H2O1coords,O1coords,O2coords,q.xboxlength)*180/math.pi<Hbondangle  or\
                           vec.angle(H1O2coords,O2coords,O1coords,q.xboxlength)*180/math.pi<Hbondangle  or\
                           vec.angle(H2O2coords,O2coords,O1coords,q.xboxlength)*180/math.pi<Hbondangle :
                               gotHbondangle=True;
                
                startt = trackingresL[pos1][pos2][0]   
                lastt = trackingresL[pos1][pos2][1]        
                if dist1 < maxdist and gotHbondangle:                           
                    if startt == -1:                         
                         trackingresL[pos1][pos2][0] = currentt
                         trackingresL[pos1][pos2][1] = currentt
                    else:                         
                         if currentt - (lastt + q.smalleststep) < t_star: # note: dist1<maxdist already confirms there is residence in this step.
                            trackingresL[pos1][pos2][1] = currentt
                         else:#================================================================
                            if startt > q.statesL[0][0] + t_star: # remember that the first time this script actually reads is the second time
                                                                  # slot.
                                timediff = (lastt + q.smalleststep) - startt                 
                                maxbin = len(resL) - 1 
                                whichbin = timediff / q.smalleststep
                                if whichbin > maxbin:                      
                                    while maxbin < whichbin: # -1 is for the label for the pair at the first slot
                                        if len(resL) == 1:
                                            resL.append([timediff,0])
                                        else:
                                            resL.append([resL[-1][0] + q.smalleststep,0])  
                                        maxbin+=1
                                    resL[-1][1]+=1 #last entry
                                else:                       
                                    resL[whichbin][1]+=1   
                                #q.debugrestf.write("%d, %d , %d, %d ,
                                #%d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep)) #================================================================
                            trackingresL[pos1][pos2][0] = currentt
                            trackingresL[pos1][pos2][1] = currentt
                else:                     
                    if lastt != -1 and currentt - lastt >= t_star:        
                        # note: dist1>=maxdist already confirms there is no
                        # residence in this time.
#================================================================
                        if startt > q.statesL[0][0] + t_star: # remember that the first time this script actually reads is the second time
                                                              # slot.
                            timediff = (lastt + q.smalleststep) - startt                                
                            maxbin = len(resL) - 1 
                            whichbin = timediff / q.smalleststep
                            if whichbin > maxbin:                      
                                while maxbin < whichbin: # -1 is for the label for the pair at the first slot
                                    if len(resL) == 1:
                                        resL.append([timediff,0])
                                    else:
                                        resL.append([resL[-1][0] + q.smalleststep,0])  
                                    maxbin+=1
                                resL[-1][1]+=1 #last entry
                            else:                       
                                resL[whichbin][1]+=1 
                            if q.collectblockerdata:
                                blockermindist=999
                                for index in q.blockerindexL: 
                                    blockercoords=q.statesL[-1][index][q.xl:q.xl+3]
                                    blockerdist=0.5*(minimage3Ddist(blockercoords,O1coords,q.xboxlength)\
                                        +minimage3Ddist(blockercoords,O2coords,q.xboxlength))
                                    if blockerdist<blockermindist:
                                        blockermindist=blockerdist
                                whichbin = int(blockermindist/0.1)+1
                                koblib.histo(q.blockerdataL,0.1,whichbin,1)

                            if dist1 < maxdist:
                                q.lostanglecount+=1 #only angle broke
                            #q.debugrestf.write("%d, %d , %d, %d ,
                            #%d\n"%(timediff/q.smalleststep,q.statesL[-1][i][0],q.statesL[-1][j][0],trackingresL[pos1][pos2][0],trackingresL[pos1][pos2][1]+q.smalleststep)) #================================================================
                        trackingresL[pos1][pos2][0] = -1
                        trackingresL[pos1][pos2][1] = -1
                pos2+=1
            pos1+=1
      

def calcFPTRDF(q,RDFL):
    
    type1=RDFL[0][0][0]
    type2=RDFL[0][1][0]
    list1=[]
    list2=None
    currtime=q.statesL[-1][0]
    starttime=q.statesL[0][0]
    whichtimeseriesindex=-999
    whichtimebin=(currtime-starttime)/q.smalleststep-1
    count=0
    
    for item in q.FPTtypes:
        type=item[0]
        if type1==type:
            whichtimeseriesindex=count
        count+=1
    if whichtimeseriesindex==-999:
        sys.exit("can't find whichtimeseriesindex")
    
    list1=q.timeseries[whichtimeseriesindex][whichtimebin][1]
    RDFL[0][2]+=len(q.timeseries[whichtimeseriesindex][whichtimebin][1])
    
    list2 = q.Typemap[type2]
   
    for i in list1:
            #print i
            atominfo_i=q.statesL[-1][i]
            for j in list2:              
              if j != i:                       
                xdiff = minimagedist(q.statesL[-1][j][2],atominfo_i[2],q.xboxlength)
                ydiff = minimagedist(q.statesL[-1][j][3],atominfo_i[3],q.yboxlength)
                zdiff = minimagedist(q.statesL[-1][j][4],atominfo_i[4],q.zboxlength)
                dist = math.sqrt(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)
                
                              
                whichbin = int(dist / q.binsize) + 1 # +1 to account for the label slot
                koblib.histo(RDFL,q.binsize,whichbin,1)

      
def specialcompute(q):  
    if q.findFPT:    
        currtime=q.statesL[-1][0]
        q.bigFPTRDFL[0]+=1
        count = 1 #the 0th place is counter
        for pair in q.FPTpairL:                    
            calcFPTRDF(q,q.bigFPTRDFL[count])
            count+=1
    if q.calculate_residence:
        for pair in q.respairL:
            calcResT(q,q.bigtrackingresL[count],q.bigresL[count])
            count+=1
    if q.calculate_waterresidence:
        calcwaterResT(q,q.watertrackingresL,q.waterresL)
           
