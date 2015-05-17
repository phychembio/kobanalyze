from koblib import *
from kobcompute import *
from kobinfo import*

def movetominimage(r2,r1,boxwidth):
    diff = r2 - r1
    while math.fabs(diff) > boxwidth / 2:
        if diff > 0:
            r1 += boxwidth
        else:
            r1 -= boxwidth
        diff = r2 - r1
    return r1

def identifyH3O(q):
    tempOL=[]
    state=q.statesL[-1]

    OL=q.Typemap[q.OHID]+q.Typemap[q.OWID]

    for Oindex in OL:  
        tempOL.append([Oindex])

    for Hindex in q.Typemap[q.HHID]+q.Typemap[q.HWID]:                
        smallestdist=999
        minOnumber=999
        count=0
        Hcoords=state[Hindex][q.xl:q.xl+3]
        boxlength=q.xboxwidth
        for Ocount in range(len(OL)):            
            Oindex=OL[Ocount]
            Ocoords=state[Oindex][q.xl:q.xl+3]
            dist=minimage3Ddist(Hcoords,Ocoords,boxlength)
            if dist<smallestdist:
                smallestdist=dist
                minOnumber=count
            count+=1
        tempOL[minOnumber].append(Hindex)
    state[2][q.OHID]=[]
    q.Typemap[q.OHID]=[]
    count=0
    for item in tempOL:
        if len(item)==4:
            state[2][q.OHID].append(item[0])
            q.Typemap[q.OHID].append(item[0])
            break




def disttopoint(refpoint, atomindex,q,whichcec):
    atomx = q.statesL[-1][atomindex][q.xl] + q.statesL[-1][atomindex][5] *q.xboxwidth
    atomy = q.statesL[-1][atomindex][q.xl+1]+ q.statesL[-1][atomindex][6] *q.yboxwidth
    atomz = q.statesL[-1][atomindex][q.xl+2]+ q.statesL[-1][atomindex][7] *q.zboxwidth
    diffx=atomx-refpoint[0]                                                                
    diffy=atomy-refpoint[1]
    diffz=atomz-refpoint[2]
    if len(q.statesL)>1:   
        lastindex=q.statesL[-2][2]["OH"][0]        
        lastx=q.distL[whichcec][-1][1]
        lasty=q.distL[whichcec][-1][2]
        lastz=q.distL[whichcec][-1][3]
        
        diffx=movetominimage(lastx,diffx,q.xboxwidth) #move to the same cell to where last x,y,z were
        diffy=movetominimage(lasty,diffy,q.yboxwidth)
        diffz=movetominimage(lastz,diffz,q.zboxwidth)             
        #print lastx,diffx,lasty,diffy,lastz,diffz    
    dist=math.sqrt(diffx**2+diffy**2+diffz**2)   
    return (diffx,diffy,diffz, dist)

def calcdistances(q):
    lastH3OID=None
       
    for i in range(q.noofCEC):
        if len(q.statesL)>1:
            lastH3OID=q.lastH3OinfoL[i][0]-3          
        currentH3Oindex=q.statesL[-1][2].get("OH")[i]
        currentH3OID=q.statesL[-1][currentH3Oindex][0]
        currentH3Opos=q.statesL[-1][currentH3Oindex][q.xl:q.xu]
        firstH3Oindex=q.statesL[0][2].get("OH")[i]
        firstH3Opos=q.statesL[0][firstH3Oindex][q.xl:q.xu]
        tempdistL=disttopoint(firstH3Opos,q.Typemap.get("OH")[i],q,i)        
        q.distL[i].append([currentH3OID,tempdistL[0],tempdistL[1],tempdistL[2],tempdistL[3]])
        if len(q.jumpdistL[i])==0:
            q.jumpdistL[i].append(q.distL[i][0])
        else:
            if lastH3OID!=currentH3OID: #reacted
                q.jumpdistL[i].append(q.distL[i][-1])
            else:
                q.jumpdistL[i].append([currentH3OID]+q.jumpdistL[i][-1][1:])
        orderedIDL=q.orderedIDL[currentH3OID][1]
        H3Oinfo=[orderedIDL]+q.statesL[-1][currentH3Oindex]
        q.H3OinfoL[i].append(H3Oinfo)
        if len(q.statesL)>1:               
            lastH3Opos=q.lastH3OinfoL[i][1:] #needs to modify for multiple CEC
            q.OOdistL[i].append(AIMDlib.minimage3Ddist(currentH3Opos,lastH3Opos,q.boxlength))
            #q.OOdistL.append([q.time/1000]+[minimage3Ddist([0,0,0],lastH3Opos,q.boxlength)])        
        else:
            q.OOdistL[i].append(0)

def printresults(prefix,q):
    outf=open(prefix+"distinfo.csv","w")
    for item in q.distL:
        outf.write("%f , %d , %f , %f , %f , %f\n"%(item[0],item[1],item[2],item[3],item[4], item[5]))
    outf.close()

    outf=open(prefix+"jumpdistinfo.csv","w")
    for item in q.jumpdistL:
        outf.write("%f , %d , %f , %f , %f , %f\n"%(item[0],item[1],item[2],item[3],item[4], item[5]))
    outf.close()
    outf=open(prefix+"OOdistinfo.csv","w")
    for item in q.OOdistL:
        outf.write("%f , %f \n"%(item[0],item[1]))
    outf.close()
    outf=open(prefix+"H3Oindexinfo.csv","w")
    for item in q.H3OindexL:
        outf.write("%f , %d, %f, %f, %f\n"%(item[0],item[1], item[2], item[3],item[4]))
    
    #outf=open(prefix+"S1.csv","w")
    #for item in q.S1histL[1:]:
    #    outf.write("%d, %d\n"%(item[0],item[1]))
    #outf.close()

    #outf=open(prefix+"S2.csv","w")
    #for item in q.S2histL[1:]:
    #    outf.write("%d, %d\n"%(item[0],item[1]))
    #outf.close()

    #outf=open(prefix+"ringsizes.csv","w")
    #for item in q.ringsizehistL[1:]:
    #    outf.write("%d, %d\n"%(item[0],item[1]))
    #outf.close()

    #outf=open(prefix+"noofrings.csv","w")
    #for item in q.numberinghistL[1:]:
    #    outf.write("%d, %d\n"%(item[0],item[1]))
    #outf.close()


    outf.close()