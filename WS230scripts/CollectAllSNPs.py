# script to collect data into single file

import sys
#import re

def addCompositeRecord(fullpath,basedir,expname):
	f=open(fullpath,'r')
	outfile=basedir+"/vars_"+expname+".txt"
	o=open(outfile,'w')
	lines=f.readlines()
	o.write("\t".join(["ID","Chr","Pos","N2var","CBvar",expname+"_N2fwd",expname+"_N2rev",expname+"_CBfwd",expname+"_CBrev",expname+"_N2counts",expname+"_CBcounts",expname+"_readDepth",expname+"_CBfreq","\n"]))	
	for line in lines:
		if line[0]=="#":
			continue
		data=line.split("\t")
		id=data[0]+"_"+data[1]
		counts=[]		
		if len(data)>8:
			info=data[7]
			I16=info.split(";")[1]
			if "I16=" in I16:
				counts=map(lambda x: int(x), I16.split("=")[1].split(",")[0:4])
			else:
				print("No I16 tag")
				counts=[0,0,0,0]
		else:
			print("No Info column")
			counts=[0,0,0,0]
		N2_counts=counts[0]+counts[1]
		CB_counts=counts[2]+counts[3]
		read_depth=sum(counts)
		if read_depth!=0:
			CB_freq=float(CB_counts)/read_depth
		else:
			CB_freq='NA'
		counts=counts+[N2_counts,CB_counts,read_depth,CB_freq]
		data.insert(0,id)
		newdata=data[0:3]+data[4:6]+map(lambda x: str(x),counts)
		o.write("\t".join(newdata)+"\n")		
	f.close()
	o.close()

# def addCompositeName(fullpath,basedir,expname):
# 	f=open(fullpath,'r')
# 	expname=filename.split("a.")[0]
# 	outfile=expname+"b.txt"
# 	o=open(outfile,'w')
# 	lines=f.readlines()
# 	for line in lines:
# 		data=line.split("\t")
# 		id="CHROMOSOME_"+data[3]+"_"+data[4]
# 		data.insert(0,id)
# 		o.write("\t".join(data))		
# 	f.close()
# 	o.close()

if __name__=='__main__':
	for fullpath in sys.argv[1:]:
		print("processing "+fullpath)
		i=fullpath.rfind('/')
		filename=fullpath[i+1:]
		expname='_'.join(filename.split('_')[1:3])
		basedir=fullpath[:i]
		basedir=basedir.replace('vcfFiles','finalData')
		addCompositeRecord(fullpath,basedir,expname)

