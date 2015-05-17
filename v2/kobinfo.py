#!/usr/bin/python

class KobInfo:
    
    stepskipped=0
    needCEC=False
    MSDtypes,pairL,respairL,closestL,notpolymerL, atommapC1S,OH_L,lastorderL,orderL,counthopsL,bigMSDL,bigRDFL,bigClosestL,bigresL,hopsL,C1SL,statesL,waterorderL,molL, bigciL,NstateL = ([] for i in range(21))
    Typemap=[[]] #first list is empty so that the the indexes match the actual atom types
    def __init__(self,trajfile,evbfile,trajext,evbext,molfile,restartfile,smalleststep,binsize,noofCEC,skipstep,needmolID,computeC1Sdist,restart,MSDdecomposition,calculateclosest,c_i_info):
        self.trajfile=trajfile
        self.evbfile=evbfile
        self.trajext=trajext
        self.evbext=evbext
        self.molfile=molfile
        self.restartfile=restartfile
        self.calculateclosest=calculateclosest        
        self.smalleststep=smalleststep
        self.binsize=binsize
        self.noofCEC=noofCEC # no need to specify if needCEC is true (enabled by having -1 in the above lists)
        self.skipstep=skipstep #time steps need to skip 
        self.needmolID=needmolID  #needed for reactive dynamics        
        self.computeC1Sdist=computeC1Sdist
        self.restart=restart        
        self.MSDdecomposition=MSDdecomposition
        self.c_i_info=c_i_info
        

