import math

class KobInfo:
    
    stepskipped=0
    needCEC=False
    noofCEC=-1
    smalleststep=0
    noofatom=0
    nooftype=0        
    anglebinsize=0
    interface=False
    needmolID=False 
    correctdrift=False#needed for reactive dynamics
    computeC1Sdist=False
    MSDdecomposition=False
    calculateclosest=False
    calculate1DSqdisp=False
    c_i_info=False
    calculate_residence=False
    calculate_waterresidence=False
    calculate_endtoenddist=False
    findwaterclusters=False
    findFPT=False
    calculateFPTproperties=False
    binsize=0.05
    ci_binsize=0.025
    anglebinsize=1.25*(4*math.atan(1))/180
    skipstep=0
    trajext='.lammpstrj'
    evbext='.evb'

    
    MSDtypes,FPTpairL,bigFPTL,bigtrackingFPTL, angletypes,pairL,CECorderL,bigMSDL,bigangleL,bigRDFL,watersizeclustersL,nowaterclustersL,statesL, Typemap = ([] for i in range(14))
 
# essential lists 
   #statesL
   #MSDtypes
   #pairL
   #CECorderL
   #bigMSDL
   #bigRDFL
   #Typemap
#not essential lists  
    #respairL,closestL,conditionL,notpolymerL, atommapC1S,OH_L,counthopsL,
    #bigClosestL,bigConditionalL,bigtrackingresL,bigresL,hopsL,C1SL,statesL,waterorderL,molL, bigciL,NstateL