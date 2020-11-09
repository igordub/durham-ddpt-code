#!/usr/bin/python

#Workflow based on Trajectory Converter - v1.5 by: Justin Lemkul 
#completely reimplemented and improved by Peter Schmidtke & Jesus Seco
#updated by Tom Rodgers

import sys,os,re,fnmatch

if len(sys.argv)>4 :
	f=sys.argv[1]
	if not os.path.exists(f):
		sys.exit(" ERROR : Something you provided does not exist. Breaking up.\n\nUSAGE : python trajconv_peter.py amberCrd amberTop trajDir trajPattern outPutPrefix\n\n \
Example : python amber2xtc.py mdcrd.crd mdcrd.top md *.x.gz md_gromacs\n")
else :
	sys.exit(" \n  USAGE : python amber2xtc.py AMBERCRD AMBERTOP TRAJDIR TRAJPATTERN OUTPUTPREFIX\n\
      Example : python amber2xtc.py mdcrd.crd mdcrd.top md *.x.gz md_gromacs\n\
      Note that the AmberCrd can also be a PDB file.\n")

crd=sys.argv[1]
top=sys.argv[2]
trajdir=sys.argv[3]
pattern=sys.argv[4]
outputPref=sys.argv[5]

traj_files=fnmatch.filter(os.listdir(trajdir),pattern) #get the fpocket output folders


RE_DIGIT = re.compile(r'(\d+)')     #set a pattern to find digits
ALPHANUM_KEY = lambda s: [int(g) if g.isdigit() else g for g in RE_DIGIT.split(s)]      #create on the fly function (lambda) to return numbers in filename strings
traj_files.sort(key=ALPHANUM_KEY)      #sort by these numbers in filenames

print "Will convert the following files : "
print traj_files

csn=1
for file in traj_files :
	ptrajtmp=open("ptraj_tmp.ptr","w")
	print "currently converting "+file
	ptrajtmp.write("trajin "+trajdir+os.sep+file+"\n")
	ptrajtmp.write("reference "+crd+"\n")
	ptrajtmp.write("center :120-130 mass origin\n")
	ptrajtmp.write("image origin center familiar\n")#:* byres
	
	ptrajtmp.write("trajout pdb_tmp/mdcrd.pdb pdb\n")
	ptrajtmp.write("go")	
	ptrajtmp.close()
	if not os.path.exists("pdb_tmp"):
		os.mkdir("pdb_tmp")

	os.system("ptraj "+top +"< ptraj_tmp.ptr")

	if not os.path.exists("xtc_tmp"):
		os.mkdir("xtc_tmp")
	#move to *.pdb
	os.system("cd pdb_tmp; for i in mdcrd.pdb.*; do a=${i##*.}; mv mdcrd.pdb.${a} mdcrd_${a}.pdb; done ; cd ..")
	
	pdb_files=fnmatch.filter(os.listdir("pdb_tmp"),"*.pdb")
	pdb_files.sort(key=ALPHANUM_KEY)      #sort by these numbers in filenames	
	if csn==1:
		os.system("editconf_d -f pdb_tmp/mdcrd_1.pdb -o "+outputPref+"_t1_top.gro")

	for pdb in pdb_files:
		os.system("echo \"0\" | trjconv_d -s pdb_tmp/"+pdb+" -f pdb_tmp/"+pdb+" -o xtc_tmp/traj_"+str(csn)+".pdb.xtc -t0 "+str(csn))
		csn+=1

	if os.path.exists(outputPref+"_traj.xtc"):
		os.system("trjcat_d -f "+outputPref+"_traj.xtc xtc_tmp/*.pdb.xtc -o "+outputPref+"_traj.xtc")
	else :
		os.system("trjcat_d -f xtc_tmp/*.pdb.xtc -o "+outputPref+"_traj.xtc")

	os.system("rm -rf pdb_tmp/*.pdb")	
	os.system("rm -rf xtc_tmp/*.xtc")


os.remove("ptraj_tmp.ptr")
os.system("rmdir pdb_tmp")
os.system("rmdir xtc_tmp")

