import os,sys

# FILE LISTS COME FROM PUBS DB TOOLS
# specifically, pubs/dlleepubs/utils/dump_
# ------------------------------------------------------------------------

flists = ["training_file_list.txt","testing_file_list.txt"]
#flists = ["filelists/cocktail_mcc8v4/dblist_mcc8v4_cocktail_p00.txt",
#          "filelists/cocktail_mcc8v4/dblist_mcc8v4_cocktail_p01.txt",
#          "filelists/cocktail_mcc8v4/dblist_mcc8v4_cocktail_p02.txt",
#          "filelists/cocktail_mcc8v4/dblist_mcc8v4_cocktail_p03.txt",
#          "filelists/cocktail_mcc8v4/dblist_mcc8v4_cocktail_p04.txt"]
os.system("mkdir inputlists")

joblist = open("joblist.txt",'w')

njobs = 0
for f in flists:
    fin = open(f,'r')
    fls = fin.readlines()
    for l in fls:
        l = l.strip()
        i = l.split('_')
        run = str(i[0])
        subrun = str(i[1])
        #jobid = 10000*run + subrun
        jobid =  subrun
        print >> joblist,jobid
        infile = open('inputlists/input_%s.txt'%(jobid),'w')
        print >> infile,'.'.join([str(l),'root'])
        #print >> infile
        infile.close()
        njobs+=1
print "number of jobs: ",njobs
joblist.close()
        
        


