import os,sys
import ROOT as rt

outdir="/cluster/kappa/90-days-archive/wongjiradlab/twongj01/ssnet_retrain_cropper/output/cocktail/p03"

# read in joblist
#fjob = open("joblist.txt",'r')
fjob = open("rerunlist.txt",'r')
ljob = fjob.readlines()
fjob.close()

jobids = []
for l in ljob:
    l = l.strip()
    jobids.append( int(l) )

print "num jobids: ",len(jobids)


completelist = []
missinglist = []
for jobid in jobids:
    rfile = "%s/ssnet_retraining_%06d.root"%(outdir,jobid)
    if not os.path.exists(rfile):
        missinglist.append(jobid)
        continue
    ttree = rt.TChain("image2d_adc_tree")
    ttree.Add(rfile)
    nentries = 0
    try:
        nentries = ttree.GetEntries()
    except:
        missinglist.append(jobid)
        continue

    if nentries>0:
        completelist.append(jobid)
    else:
        missinglist.append(jobid)

print "jobs to still run: %d/%d"%(len(missinglist),len(jobids))

frerun = open("rerunlist.txt", 'w')
for jobid in missinglist:
    print >> frerun,jobid
frerun.close()
        


