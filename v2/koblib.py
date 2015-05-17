#!/usr/bin/python
import re
import sys
import math
import kobcompute

def initialize(q):
    
    if q.computeC1Sdist:
       q.C1SL=[0]


    for type in q.MSDtypes:
        q.bigMSDL.append([type])
        if type==-1:
            q.needCEC=True
        if type==-2 or type==-3:
            q.MSDdecomposition=True

    q.bigRDFL.append(0)
    for pair in q.pairL:
        q.bigRDFL.append([pair])
        if -1 in pair[0] or -1 in pair[1]:
            q.needCEC=True

    for pair in q.respairL:
        q.bigresL.append([(pair[0],pair[1])]) 
        q.bigresL[-1].append([0,0])#[0,0]=[timediff=0,count=0]
        if -1 in pair[0] or -1 in pair[1]:
            q.needCEC=True

    q.bigClosestL.append(0)

    for pair in q.closestL:
        if -1 in pair[0] or -1 in pair[1]:
            q.needCEC=True

    for item in q.closestL:
        q.bigClosestL.append([item])

    if  q.MSDdecomposition:
        if len(hopsL)==0:
           for i in range(0,q.noofCEC):
              q.hopsL.append([])
              q.counthopsL.append([0,0,[],[]])
    if q.c_i_info:
        q.bigciL.append(0)        
        q.bigciL.append([1])
        q.bigciL.append([2])
        q.NstateL.append(["Nstate"])

def findprintfreq(trajfile):
    trajf=open(trajfile,"r")
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time1=int(trajf.next())  
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time2=int(trajf.next())
    trajf.close()
    return time2-time1



def histo(binL,binsize,whichbin):
     maxbin=len(binL)-1 
     if whichbin>maxbin:                      
            while maxbin<whichbin: # -1 is for the label for the pair at the first slot                          
                if len(binL)==1:
                  binL.append([binsize/2,0])
                else:
                  binL.append([binL[-1][0]+binsize,0])  
                maxbin+=1
            binL[-1][1]+=1 #last entry
     else:
            binL[whichbin][1]+=1    

def setupTypemap(statesL,Typemap):
    Typemap.append([])
    index=3
    currentatomType=1
    for item in statesL[0][3:]:    
        if item[1]==-1:
            print "here"
        if item[1]>(len(Typemap)-1):        
            diff=(item[1]+1)-len(Typemap)
            for i in range(0,diff):
                Typemap.append([])
            currentatomType=item[1]  
        elif item[1]==-1 and currentatomType!=-1:  
            Typemap.append([])
            currentatomType=item[1]  
        Typemap[currentatomType].append(index)
        index+=1

def readmolID(molfile):        
        molf=open(molfile,'rU')
        for i in range(0,9):
            line=molf.next()
        for line in molf:
            info=line.split()
            atomid=int(info[0])
            typeid=int(info[1])
            molid=int(info[2])        
            if typeid==OWID or typeid==OHID:
                q.molL.append([atomid,molid])
        molf.close()
        return sorted(q.molL, key=lambda k:k[1])

        if q.needmolID: q.molL=readmolID() 


def gotoline(fstream, text, after):
    
    for line in fstream:
        match=re.match(text,line)
        if match:
           found= True
           while after > 0:        
              line=fstream.next()
              after=after -1
           return line
    return 'String Not Found'   



def inputstate(trajf,timestep):
    state=[]  #(timestep, noofatom, halfwidth, (atom_1)...(atom_N)(CEC_1)..(CEC_N))
    #atom/CEC info=(atomid, atom type, x, y, z)
    
    state.append(timestep)
    gotoline(trajf,'ITEM: NUMBER',0)
    noofatom=int(trajf.next())
    state.append(noofatom)
    gotoline(trajf,'ITEM: BOX',0)
    xinfo=trajf.next().split()
    yinfo=trajf.next().split()
    zinfo=trajf.next().split()
    
    xlo=float(xinfo[0])
    xhi=float(xinfo[1])
    ylo=float(yinfo[0])
    yhi=float(yinfo[1])
    zlo=float(zinfo[0])
    zhi=float(zinfo[1])
    boxwidth=(float(xinfo[1])-float(xinfo[0]))
    state.append(boxwidth/2)
    line=gotoline(trajf,'ITEM: ATOMS',0)
    xorxs=line.split()[5]
    

    tempL=[]
    for i in range(0,noofatom):
        atomdata=trajf.next().split()    
        atomtype=atomdata[1] 
        if atomtype!="98":
            if xorxs=="x":
                x=float(atomdata[3])
                y=float(atomdata[4])
                z=float(atomdata[5])
            else:
                x=float(atomdata[3])*boxwidth-xlo
                y=float(atomdata[4])*boxwidth-ylo
                z=float(atomdata[5])*boxwidth-zlo          
            tempL.append((int(atomdata[0]),int(atomdata[1]),x,y,z))
    for i in sorted(tempL,key=lambda k: (k[1],k[0])):
        state.append(i)
    return state

def fixpbc(tempstate, statesL):
      for i in range(3,len(tempstate)):
               if len(statesL)==0:
                  tempstate[i]=tempstate[i]+(0,0,0) 
               else:    
                        
                  crossxdist=tempstate[i][2]-statesL[-1][i][2]
                  crossydist=tempstate[i][3]-statesL[-1][i][3]
                  crosszdist=tempstate[i][4]-statesL[-1][i][4]
                  crossx=statesL[-1][i][5]
                  crossy=statesL[-1][i][6]                  
                  crossz=statesL[-1][i][7]
                  hboxwidth=tempstate[2]
                  if abs(crossxdist) >= hboxwidth:
                     if crossxdist< 0:
                        crossx = crossx + 1
                     else:
                        crossx = crossx - 1                  
                  if abs(crossydist) >= hboxwidth:
                     if crossydist< 0:
                        crossy = crossy + 1
                     else:
                        crossy = crossy - 1
                  if abs(crosszdist) >= hboxwidth:
                     if crosszdist< 0:
                        crossz = crossz + 1
                     else:
                        crossz = crossz - 1            		                                     		                                     		                                     		             
                  tempstate[i]=tempstate[i]+(crossx,crossy,crossz)     

                  
                
def insertCEC( state, evbf, waterorderL,nextorderL,CECtimestep,q):          
          
          pivotL=[] #stores possible new mol is for CEC  
          whichCEC=-999          
          q.OH_L=[]
          CEC_coordsL=[]
          discretCEC_coordsL=[]
          contCEC_coordsL=[]
          firststep=False
          q.noofCEC=int(evbf.next().split()[1])
          evbf.next()
          if len(nextorderL)==0:
            firststep=True
            for i in range(0, q.noofCEC):                                  
                 mol_id=int(evbf.next().split()[1])                 
                 nextorderL.append([mol_id,i,0]) 
                 q.OH_L.append(mol_id)
 #---------------------------           
            if q.needmolID:
               tempOW_L=[]
               for item in q.molL:
                   if item[1] not in q.OH_L:
                      tempOW_L.append(item)
               count=1
               for item in sorted(tempOW_L,key=lambda k: k[0]):
                   nextwaterorderL.append([item[0],item[1],count])
                   count+=1
 #-----------------------------          
   


  
          currentorderL=[]          
          for s in nextorderL:
            currentorderL.append(s[:])
 #---------------------------           
          if q.needmolID:
              currentwaterorderL=[]
              for s in nextwaterorderL:
                  currentwaterorderL.append(s)
 #-----------------------------        

          if CECtimestep==state[0]:
             for item in nextorderL:
                item[2]=0   
       
         
          for line in evbf:                      
             match1=re.search('COMPLEX (\d+): (\d+)',line)
             match2=re.search('NEXT_PIVOT_STATE ([0-9]\d*)',line)
             match3=re.search('CEC_COORDINATE',line)
             match4=re.search('END_OF_COMPLEX ' +str(q.noofCEC),line)                   
             
             if match1:
                pivotL=[]
                whichCEC=int(match1.group(1))-1
                noofstates=int(match1.group(2))
                evbf.next()                       
                for i in range(0,noofstates):
                   pivotL.append(int(evbf.next().split()[4]))            
             if match2:                
                pivotL_index=int(match2.group(1))
                
                if pivotL_index!=0:
#------------------------------------           
                    if q.needmolID: #swap the current hydronium to water for next step                        
                        for i in range(0,len(q.nextwaterorderL)):
                           if q.nextwaterorderL[i][1]==pivotL[pivotL_index]:
                                # print("%d -->%d"%(q.nextwaterorderL[i][1],q.nextorderL[whichCEC][0]))
                                 q.nextwaterorderL[i][1]= q.nextorderL[whichCEC][0]                                 
                                 q.nextwaterorderL[i][0]= q.molL[q.nextwaterorderL[i][1]-q.molL[0][1]][0]
                                 break                         
#--------------------------------
                    nextorderL[whichCEC][0]=pivotL[pivotL_index] #change mol_id for next CEC      
                    nextorderL[whichCEC][2]=1 #reacted
          

             if match3:                
                coordinates=[float(x) for x in evbf.next().split()]
               
                if CECtimestep==state[0]:                   
                   CEC_coordsL.append((currentorderL[whichCEC][1],currentorderL[whichCEC][0],-1,coordinates[0],coordinates[1],coordinates[2]))

             if q.c_i_info:
                    match5=re.search('EIGEN_VECTOR',line)          
                    if match5:
                        line=evbf.next();                        
                        kobcompute.collect_ci(sorted([float(x) for x in line.split()],reverse=True)[:2],q,noofstates)
                 
                
             if match4:                         
                 if CECtimestep==state[0]:
                    for coords in sorted(CEC_coordsL):
                        state.append(coords[1:])      
#--------------------------------
                    if q.needmolID:                    
                        mapkeyL=[]
                        firstwaterindex=3
                        for atominfo in state[3:]:
                            if atominfo[1]==OWID:    
                                break
                            firstwaterindex+=1                        
                        for item in currentwaterorderL:
                           index=firstwaterindex
                           for atominfo in state[firstwaterindex:]:
                               if item[0]==atominfo[0]:
                                   mapkeyL.append([item[2],index])
                                   break
                               index+=1
                        waterL=[]    
                        for item in sorted(mapkeyL, key=lambda k: k[0]):
                            waterL.append(state[item[1]])
                        count=0
                        for item in sorted(mapkeyL, key=lambda k: k[1]):
                            state[item[1]]=waterL[count]
                            count+=1
                 return (sorted(nextorderL),sorted(currentorderL,key=lambda k: k[1]))                 
#--------------------------------                 
          
            
def evbtime(evbf):
     line=gotoline(evbf, 'TIMESTEP ',0)
     if line!="String Not Found":     
        return int(line.split()[1])
     else:
        print "End of file: evbf"
        while True:
            evbf.next()
        

def trajtime(trajf,q):                       
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time=int(trajf.next()) 
    
    while q.stepskipped<q.skipstep:
        line=gotoline(trajf,'ITEM: TIMESTEP\n',0);
        time=int(trajf.next()) ;   
        q.stepskipped+=1

    if line!="String Not Found":     
        return time 
    else:
        print "End of file: trajf"
        while True:
            trajf.next()
        

def eventime(evbf,trajf,q):
    ttime=trajtime(trajf,q)
    etime=evbtime(evbf)
    while etime>ttime:
        ttime=trajtime(trajf,q)
    while ttime>etime:
        etime=evbtime(evbf)
    return ttime

def countatoms(state, atomIDlist):
    count=0
    for atom in state[3:]:
        if atom[1] in atomIDlist:
            count+=1
    return count

def dotvec(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def lenvec(v):
    return math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])



def printresults(args,q):
    boxwidth=q.statesL[-1][2]*2
    for i in range(0,len(q.MSDtypes)):
        msdf=open('MSD_type='+str(q.MSDtypes[i])+'_'+args+'.txt','w')
        for MSD in q.bigMSDL[i][1:]:
            msdf.write("%d   %f \n" %( MSD[0] , MSD[1]/MSD[2]))
        msdf.close()

    for i in range(0,len(q.respairL)):
        resf=open('ResT_type='+str(q.respairL[i][0])+'-'+str(q.respairL[i][1])+'-'+str(q.respairL[i][2])+'_'+args+'.txt','w')
        for resL in q.bigresL[i][1:]:
            resf.write("%d   %f \n" %( resL[0] , resL[1]))
        resf.close()


    count=1  #the 0th place is counter
    
    for pair in q.pairL:        
        nooftype1=countatoms(q.statesL[-1],pair[0])
        nooftype2=countatoms(q.statesL[-1],pair[1])
        RDF=q.bigRDFL[count]
        rdff=open('RDF_type='+str(pair[0])+'_'+str(pair[1])+args+'.txt','w')
        sum=0
        for bin in RDF[1:]:
            sum+=bin[1]            
            rdff.write("%f   %f  %f\n" %(bin[0] , bin[1]/(4*math.pi*bin[0]**2*q.binsize*nooftype2/boxwidth**3*nooftype1*q.bigRDFL[0]),sum*1./q.bigRDFL[0]/nooftype1))
        rdff.close() 
        count+=1
    
    if q.calculateclosest==True:
        count=1
        for pair in q.closestL:            
            nooftypeatoms=countatoms(q.statesL[-1],pair[0])
            if pair[1]==999:
                towhere="toPolymer"
            else:
                towhere=str(pair[1])
                
            closestf=open("closest_type="+str(pair[0])+"_"+towhere+args+'.txt','w')
            subcloseL=q.bigClosestL[count]
            for bin in subcloseL[1:]:
                closestf.write("%f   %f \n" %(bin[0] , bin[1]*1./(q.bigClosestL[0]*nooftypeatoms)))
            closestf.closse()
            count+=1
    if q.computeC1Sdist==True:
        count=1
        nooftypeatoms=countatoms(q.statesL[-1],[C1ID])     
                
        closestf=open("C1S"+args+'.txt','w')
        
        for bin in q.C1SL[1:]:
            closestf.write("%f   "%bin[0])  
            for item in bin[1:]:
               closestf.write("%d  " %(item))
            closestf.write("\n")
        closestf.close()
        count+=1
    if q.MSDdecomposition:
        hopf=open("counthops.txt","w")
        hopMSDf=open("hopMSDs.txt","w")
        count=1
        for item in q.hopsL:
            print str(count) +" ",
            print item
            count+=1
        count=1
        for item in q.counthopsL:
            hopf.write("%d %d %d\n"%(count,item[0],item[1]))
            for list in item[2:]:
                for listitem in list:
                    hopMSDf.write(str(listitem)+" ")
                hopMSDf.write("\n")
            count+=1
        hopf.close()
        hopMSDf.close()
    if q.c_i_info:
        for i in range(1,3):
            L=q.bigciL[i]
            fout=open('c'+str(i)+'_'+args+'.txt','w')            
            for bin in L[2:]:
                 fout.write("%f   %f  %f\n" %(bin[0] , bin[0]*bin[0],bin[1]*1./q.bigciL[0]))
            fout.close() 

        L=q.NstateL
        fout=open('Nstate_'+args+'.txt','w')            
        for bin in L[2:]:
                fout.write("%d   %f\n" %(bin[0] ,bin[1]*1./q.bigciL[0]))
        fout.close() 