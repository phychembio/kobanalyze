import sys
import re
inf=open("out.lammpstrj","r")
outf=open("start.lammpstrj","w")
line=inf.next()
outf.write(line)
for line in inf:
    matchend=re.search("ITEM: TIMESTEP",line)
    if matchend:
        sys.exit()
    info=line.split()
    if len(info)==9:
        outf.write("%s %s %s\n"%(info[0],info[1],info[2]))
    else:
	    outf.write(line)
	
    
