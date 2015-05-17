import re
import sys
import math
import kobcompute

def initialize(q):
    gettrajinfo(q)


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


    for type in q.angletypes:     
        tempL=[]   
        for i in q.angletimeL:
            tempL.append([type])
        q.bigangleL.append(tempL)

   
        

    if q.calculateclosest:
        for pair in q.closestL:
            if -1 in pair[0] or -1 in pair[1]:
                q.needCEC=True
  
    if q.needCEC:
        getevbinfo(q)
        q.OorderL=[]
        q.OH_L=[]       
       
    
    if q.calculateclosest:
        q.bigClosestL=[0]
        q.bigConditionalL=[0]                
        count=0
        for item in q.closestL:
            q.bigClosestL.append([item])
            if item[2]==True:
                q.bigConditionalL.append([q.conditionL[count]])
            count+=1

    if  q.MSDdecomposition:
        q.hopsL=[]
        for i in range(0,q.noofCEC):
           q.hopsL.append([])
           q.counthopsL.append([0,0,[],[]])
   
             
             
    if q.c_i_info:        
        q.bigciL=[0]
        q.NstateL=([["Nstate"]])
        q.bigciL.append([1])
        q.bigciL.append([2])
        q.bigciL.append([3])
    
    if q.calculate1DSqdisp:
        q.big1DSqdisp=[0]    
        for type in q.Sqdisp:
            q.big1DSqdisp.append([type])  
    
 
    if q.findwaterclusters:
        q.nowaterclustersL=["Label"] #not utilized yet
        q.watersizeclustersL=["Label"] #not utilized yet



    if q.calc_reactivewater: q.moltoIDL,q.IDtomolL =readmolID(q) 

    


def gettrajinfo(q):
    trajfile=q.trajfile+q.trajext
    trajf=open(trajfile,"r")
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time1=int(trajf.next())  
    
    trajf.next()
    
    q.noofatom=int(trajf.next())
    
    trajf.next()
    info=trajf.next().split()
    q.xlo=float(info[0])
    q.xhi=float(info[1])
    info=trajf.next().split()
    q.ylo=float(info[0])
    q.yhi=float(info[1])
    info=trajf.next().split()
    q.zlo=float(info[0])
    q.zhi=float(info[1])
    q.boxlengthL=[]
    q.boxlengthL.append(q.xhi-q.xlo)    
    q.boxlengthL.append(q.yhi-q.ylo)
    q.boxlengthL.append(q.zhi-q.zlo)
    q.xboxlength=q.boxlengthL[0]
    q.yboxlength=q.boxlengthL[1]
    q.zboxlength=q.boxlengthL[2]
    q.boxlength=q.xhi-q.xlo
    labelL=trajf.next().split()[2:] #ITEM: ATOMS id type....
    typelabelindex=None
    for i in range(len(labelL)):
        label=labelL[i]
        if label=="type":
            typelabelindex=i
            break
    maxtypeid=0
    for line in range(q.noofatom):
        line=trajf.next()
        typeid=int(line.split()[typelabelindex])
        if typeid>maxtypeid:
            maxtypeid=typeid    
    q.nooftype=maxtypeid
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time2=int(trajf.next())
    trajf.close()
    
    q.smalleststep=time2-time1
    q.angletimeL=[q.smalleststep,5*q.smalleststep,10*q.smalleststep,50*q.smalleststep,100*q.smalleststep,500*q.smalleststep,1000*q.smalleststep,5000*q.smalleststep,10000*q.smalleststep]
    
def getevbinfo(q):
    evbfile=q.evbfile+q.evbext
    evbf=open(evbfile,"r")
    line=gotoline(evbf, 'TIMESTEP ',0)
    q.noofCEC=int(evbf.next().split()[1])

   
def jifhisto(binL,binsize,whichbin,increment):
     maxbin=len(binL)-1 
     finalbin=whichbin
     if whichbin>maxbin:                      
        while maxbin<whichbin: 
            if len(binL)==1:
                binL.append([binsize/2,0,0])
            else:
                binL.append([binL[-1][0]+binsize,0,0])
            maxbin+=1
        finalbin=-1
     binL[finalbin][1]+=increment #last entry
     binL[finalbin][2]+=1
     

def histo(binL,binsize,whichbin,increment):
     maxbin=len(binL)-1 
     finalbin=whichbin
     if whichbin>maxbin:                      
        while maxbin<whichbin: 
            if len(binL)==1:
                binL.append([binsize/2,0,0])
            else:
                binL.append([binL[-1][0]+binsize,0,0])
            maxbin+=1
        finalbin=-1
     binL[finalbin][1]+=increment #last entry
     binL[finalbin][2]+=1
     
def setupTypemap(state,q):    
    q.Typemap=[]
    Typemap=q.Typemap
    for i in range(q.nooftype+1): #+1 is for the spaceholder to match the atomtype IDs
        Typemap.append([])
    if q.needCEC:
        Typemap.append([])
        if q.MSDdecomposition:
            Typemap.append([])
            Typemap.append([])
    index=3
    for item in state[3:]:    
        typeID=item[1]
        atomID=item[0]
        Typemap[typeID].append(index)        
        index+=1

    if q.needCEC:
        for i in range(q.noofCEC):
            if q.MSDdecomposition:               
               Typemap[-2].append(q.noofatom+q.noofCEC+i+3)
               Typemap[-3].append(q.noofatom+2*q.noofCEC+i+3)            
            Typemap[-1].append(q.noofatom+i+3)

            #   Typemap[-3].append(q.noofatom+i+3)
            #   Typemap[-2].append(q.noofatom+q.noofCEC+i+3)
            #   Typemap[-1].append(q.noofatom+2*q.noofCEC+i+3)
            #else:
            #   Typemap[-3].append(q.noofatom+i+3)


        
        

def readmolID(q):        
        
        molf=open(q.molfile,'rU')
        for i in range(0,9):
            line=molf.next()
        for line in molf:
            info=line.split()
            atomid=int(info[0])
            typeid=int(info[1])
            molid=int(info[2])        
            if typeid==q.OWID or typeid==q.OHID:
                q.moltoIDL.append([molid,atomid])
                q.IDtomolL.append([atomid,molid])
        molf.close()
        return [sorted(q.moltoIDL),sorted(q.IDtomolL)]

        


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



def inputstate(trajf,q,timestep):
    q.typeL=[-1]    
    q.coordsL=[[-999,-999,-999]]
    state=[]  #(timestep, q.noofatom, halfwidth, (atom_1)...(atom_N)(CEC_1)..(CEC_N))
    #atom/CEC info=(atomid, atom type, x, y, z)
    
    state.append(timestep)
    gotoline(trajf,'ITEM: NUMBER',0)
    
    state.append(q.noofatom)
    gotoline(trajf,'ITEM: BOX',0)
    
    state.append("placeholder") #will be used to store Typemap
    line=gotoline(trajf,'ITEM: ATOMS',0)    
    labelL=line.split()[2:]
    xorxs=None
    xindex=None
    molindex=None
    idindex=None
    typeindex=None
    for i in range(len(labelL)):
        label=labelL[i]
        if label=="x":            
            xorxs=label    
            xindex=i        
        if label=="id":
            idindex=i
        if label=="mol":
            molindex=i
        if label=="type":
            typeindex=i            
    tempL=[]
    for i in range(0,q.noofatom):
        atomdata=trajf.next().split()                       
        x=float(atomdata[xindex])
        y=float(atomdata[xindex+1])
        z=float(atomdata[xindex+2])
        atomID=int(atomdata[idindex])
        atomType=int(atomdata[typeindex])
        tempL.append((atomID,atomType,x,y,z))
        #q.coordsL.append([x,y,z])
        
   
    for i in sorted(tempL,key=lambda k:k[0]):        
        q.coordsL.append(i[2:5])
        q.typeL.append(i[1]) 
    for i in sorted(tempL,key=lambda k:k[1]):
        state.append(i)
    
    setupTypemap(state,q)   
    state[2]=q.Typemap



    return state


def fixpbc(tempstate, q):
      for i in range(3,len(tempstate)):
               if len(q.statesL)==0:
                  tempstate[i]=tempstate[i]+(0,0,0) 
               else:    
                        
                  crossxdist=tempstate[i][2]-q.statesL[-1][i][2]
                  crossydist=tempstate[i][3]-q.statesL[-1][i][3]
                  crosszdist=tempstate[i][4]-q.statesL[-1][i][4]
                  crossx=q.statesL[-1][i][5]
                  crossy=q.statesL[-1][i][6]                  
                  crossz=q.statesL[-1][i][7]                  
                  if abs(crossxdist) >= q.xboxlength/2:
                     if crossxdist< 0:
                        crossx = crossx + 1
                     else:
                        crossx = crossx - 1                  
                  if abs(crossydist) >= q.yboxlength/2:
                     if crossydist< 0:
                        crossy = crossy + 1
                     else:
                        crossy = crossy - 1
                  if abs(crosszdist) >= q.zboxlength/2:
                     if crosszdist< 0:
                        crossz = crossz + 1
                     else:
                        crossz = crossz - 1            		                                     		                                     		                                     		             
                  tempstate[i]=tempstate[i]+(crossx,crossy,crossz)     

                  
                
def insertCEC( state, evbf, OorderL,CECorderL,CECtimestep,q):          
          
          pivotL=[] #stores possible new mol is for CEC  
          whichCEC=-999                    
          CEC_coordsL=[]
          discretCEC_coordsL=[]
          contCEC_coordsL=[]
          OH_molL=[]
          evbf.next()
          evbf.next()
          if len(CECorderL)==0:          
            for i in range(0, q.noofCEC):                                  
                 mol_id=int(evbf.next().split()[1])                 
                 CECorderL.append([mol_id,i,0]) 
                 OH_molL.append(mol_id)
 #---------------------------           
            if q.calc_reactivewater:
               tempOL=[] #stores (atomid, mod ID)
               for item in q.IDtomolL:
                   tempOL.append(item)                               
               for item in tempOL:
                   OorderL.append([item[0],item[1]])  #stores (atomid, mod ID)
                   
               
 #-----------------------------          

         
          for line in evbf:                      
             match1=re.search('COMPLEX (\d+): (\d+)',line)
             match2=re.search('NEXT_PIVOT_STATE ([0-9]\d*)',line)
             match3=re.search('CEC_COORDINATE',line)
             match4=re.search('END_OF_COMPLEX ' +str(q.noofCEC),line)                   
             
             if match1:
                molBL=[]
                molAL=[]
                whichCEC=int(match1.group(1))-1
                noofstates=int(match1.group(2))
                evbf.next()                       
                for i in range(0,noofstates):
                   info=evbf.next().split()
                   molBL.append(int(info[4]))            
                   molAL.append(int(info[3]))
             if match2:                
                pivotL_index=int(match2.group(1))
                
                if pivotL_index!=0:
#------------------------------------           
                    if q.calc_reactivewater: #swap the current hydronium to water for next step     
                            
                        #findmolA
                        molALindex=molBLindex=None
                        molAmolID=molAL[pivotL_index]
                        molBmolID=molBL[pivotL_index]
                   
                        foundcount=0
                        for i in xrange(len(OorderL)):
                            if foundcount<2: 
                                molID=OorderL[i][1]
                                if molID==molAmolID: 
                                    molALindex=i;foundcount+=1
                                if molID==molBmolID: 
                                    molBLindex=i;foundcount+=1
                            else:
                                break                        
                        OorderL[molALindex],OorderL[molBLindex]=OorderL[molBLindex],OorderL[molALindex]
#--------------------------------
           
                    CECorderL[whichCEC][0]=molBL[pivotL_index] #change mol_id for next CEC      
                    CECorderL[whichCEC][2]=1 #reacted
          

             if match3:                
                coordinates=[float(x) for x in evbf.next().split()]
               
                if CECtimestep==state[0]:                   
                   CEC_coordsL.append((CECorderL[whichCEC][1],CECorderL[whichCEC][0],-1,coordinates[0],coordinates[1],coordinates[2]))

             if q.c_i_info:
                    match5=re.search('EIGEN_VECTOR',line)          
                    if match5 and CECtimestep%q.smalleststep==0:
                        line=evbf.next();                        
                        kobcompute.collect_ci(sorted([float(x) for x in line.split()],reverse=True)[:2],q,noofstates)
                 
                
             if match4:                         
                 if CECtimestep==state[0]:
                    for coords in sorted(CEC_coordsL):
                        state.append(coords[1:])      
#--------------------------------
                    if q.calc_reactivewater:#basically, the atom IDS in each type in state  are sorted
                        waterL=[]  #this is problematic when there is reactive and MSD is needed                                                                                 
                                   #this section restores the original unsorted order
                        #q.OL=q.Typemap[q.OWID]+q.Typemap[q.OHID]                        

                        for item in OorderL:                                                
                           for i in q.Typemap[q.OWID]:                               
                               if item[0]==state[i][0]: #finding where the atom is                                   
                                   waterL.append(state[i])                   
                                   break 
                               
                                   
                           


                                       
                        count=0          
                        #print CECtimestep
                        for i in q.Typemap[q.OWID]:                                                  
                            state[i]=waterL[count]                      
                            count+=1
                        
#--------------------------------                              
                 return (sorted(CECorderL)) #CECorderL is sorted based on the first entries becasue this is how the CECs are listed in the .evb file.         
                                            #when CECorderL is modified, the index whichCEC is based on this ordering.
                                            #when the CECs are actually inserted in the state, a sorting based on the second entries is used.
                                            #CECorderL[whichCEC][1] is inserted as the first entry in CEC_coordsL which is then sorted based on the first entries
    
          
            
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

def countatoms(map, atomIDlist):
    sum=0
    for atomID in atomIDlist:
        sum+=len(map[atomID])
    return sum


def dotvec(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def lenvec(v):
    return math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

def calcdrift(filename):
    infile=open(filename,"r")
    totalmass=0
    for line in infile:
       match=re.search("Masses",line)
       if match:
          break
      
    infile.next()
    massL=[]
    for line in infile:
       info=line.split()
       if len(info)!=0:
          massL.append(float(info[1]))
       else:
           break
        
    for line in infile:
       match=re.search("Atoms",line)
       if match:
          break

    tempatomL=[]
    infile.next()
    for line in infile:
       info=line.split()
       if len(info)!=0:
          tempatomL.append((int(info[0]),int(info[2])))
       else:
           break
    atomL=sorted(tempatomL)

    for line in infile:
       match=re.search("Velocities",line)
       if match:
          break

    infile.next()
    sumL=[0,0,0]
    count=0
    for line in infile:
       info=line.split()
       if len(info)==4:
           count+=1
           atomindex=int(info[0])
           atomtype=atomL[atomindex-1][1]
           atommass=massL[atomtype-1]
           totalmass+=atommass
           #print atommass
           for i in range(3):
              sumL[i]+=atommass*float(info[i+1])
       else:
           break
    for i in range(3):
        sumL[i]=sumL[i]/totalmass
    return sumL


def printresults(suffix,q):

    volume=q.xboxlength*q.yboxlength*q.zboxlength
    if (len(q.bigMSDL)>0):
        print q.bigMSDL[0]
        for i in range(0,len(q.MSDtypes)):
            msdf=open('MSD_type='+str(q.MSDtypes[i])+'_'+suffix+'.txt','w')
            for MSD in q.bigMSDL[i][1:]:
                msdf.write("%d   %f \n" %( MSD[0] , MSD[1]/MSD[2]))
            lastMSD=kobcompute.calclastMSD(q,q.MSDtypes[i])
            msdf.write("%d   %f \n" %( lastMSD[0] , lastMSD[1]))
            msdf.close()



    for listindex in range(len(q.angletypes)):
        for i in range(len(q.angletimeL)):
            anglef=open('angle_type='+str(q.angletypes[listindex])+'_t='+str(q.angletimeL[i])+"_"+suffix+'.txt','w')        
            #print q.bigMSDL[i]
            normalization=0
            if len(q.bigangleL[listindex][i])>1:
                for angle in q.bigangleL[listindex][i][1:]:                    
                    normalization+=angle[2]
                for angle in q.bigangleL[listindex][i][1:]:
                    anglef.write("%f   %f \n" %(angle[0] , angle[1]*1./normalization))       
                anglef.close()

    count=1  #the 0th place is counter    
    for pair in q.pairL:        
        nooftype1=countatoms(q.Typemap,pair[0])
        nooftype2=countatoms(q.Typemap,pair[1])
        if pair[0]==pair[1]:
            nooftype1=nooftype1-1
        RDF=q.bigRDFL[count]
        rdff=open('RDF_type='+str(pair[0])+'_'+str(pair[1])+suffix+'.txt','w')
        sum=0
        for bin in RDF[1:]:
            sum+=bin[1]
            intensity=bin[1]/(4*math.pi*bin[0]**2*q.binsize*nooftype2/volume*nooftype1*q.bigRDFL[0])
            integration=sum*1./q.bigRDFL[0]/nooftype1                        
            rdff.write("%f   %f  %f\n" %(bin[0] ,intensity,integration))
        rdff.close() 
        count+=1
    
    if q.calculateclosest==True:
        count=1
        for pair in q.closestL:            
            nooftypeatoms=countatoms(q.Typemap,pair[0])
            if pair[1]==999:
                towhere="toPolymer"
            else:
                towhere=str(pair[1])                
            closestf=open("closest_type="+str(pair[0])+"_"+towhere+suffix+'.txt','w')
            subcloseL=q.bigClosestL[count]
            for bin in subcloseL[1:]:
                closestf.write("%f   %f \n" %(bin[0] , bin[1]*1./(q.bigClosestL[0]*nooftypeatoms)))                
            closestf.close()

            if pair[2]==True:
                 type1=str(pair[1][0])+"-"+str(q.conditionL[count-1][0][0])
                 type2=str(q.conditionL[count-1][1][0])
                 conditionalf=open("condition_type="+type1+"_"+type2+suffix+'.txt','w')
                 nooftypeatoms=countatoms(q.Typemap,q.conditionL[count-1][0])
                 subcondL=q.bigConditionalL[count]
                 for bin in subcondL[1:]:
                    if bin[2]!=0:
                        conditionalf.write("%f   %f \n" %(bin[0] , bin[1]*1./bin[2]))                
                    else:
                        conditionalf.write("%f   0 \n" %bin[0])
                 conditionalf.close()
            count+=1


    if q.computeC1Sdist==True:
        count=1
        nooftypeatoms=countatoms(q.Typemap,[q.C1ID])     
                
        closestf=open("C1S"+suffix+'.csv','w')
        
        distsum=0
        angelsum=0
        for bin in q.C1SL[1:]:
            distsum+=bin[1]
            for number in bin[2:]:
                angelsum+=number            
        closestf.write(", dist , ")
        anglebin_deg=q.anglebinsize*180/(4*math.atan(1))
        for i in range(len(q.C1SL[1][2:])):
            closestf.write("%f, "%((i+1)*anglebin_deg))
        closestf.write("\n") 
        for bin in q.C1SL[1:]:                      
            closestf.write("%f,  %f,  "%(bin[0],bin[1]*1./(distsum*q.binsize)))    
            for item in bin[2:]:
               closestf.write("%f,  " %(item*1./(angelsum*anglebin_deg*q.binsize)))
            closestf.write("\n")
        closestf.close()
        count+=1

    if q.MSDdecomposition:
        hopf=open("counthops.txt","w")
        hopMSDf=open("hopMSDs.txt","w")
        #count=1
        #for item in q.hopsL:
        #    print str(count) +" ",
        #    print item
        #    count+=1
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
        for i in range(1,4):
            L=q.bigciL[i]
            fout=None
            if i!=3:
                fout=open('c'+str(i)+'_'+suffix+'.txt','w')
                normalization=0
                for bin in L[1:]:
                    normalization+=bin[1]
                normalization*=q.ci_binsize 
                for bin in L[1:]:
                    fout.write("%f   %f\n" %(bin[0],bin[1]*1./normalization))
            else:
                fout=open('ci_sqdiff'+'_'+suffix+'.txt','w')                   
                normalization=0
                for bin in L[1:]:
                    normalization+=bin[1]                          
                normalization*=q.ci_binsize
                for bin in L[1:]:                    
                    if bin[2]!=0:
                        fout.write("%f  %f\n" %(bin[0],-0.5961*math.log(bin[1]/normalization)))
            fout.close() 
        
        L=q.NstateL
        fout=open('Nstate_'+suffix+'.txt','w')            
        for bin in L[2:]:
                fout.write("%d   %f\n" %(bin[0] ,bin[1]*1./q.bigciL[0]))
        fout.close() 
    if q.calculate_endtoenddist:
        fout=open("end_to_end_dist_"+suffix+'.txt','w')
        for i in range(len(q.headtype)):
            sumavg=q.currentdist[i][0]/q.currentdist[i][2]
            sumsqavg=q.currentdist[i][1]/q.currentdist[i][2]
            fout.write("X sum %f, X^2 sum %f, count %d, <X> %f, <X^2> %f ,<X^2>-<X>^2 %f \n"\
                %(q.currentdist[i][0], q.currentdist[i][1], q.currentdist[i][2], sumavg,sumsqavg,sumsqavg-sumavg**2))
        fout.close()

