from scipy import stats
import sys


tnum=0
all={}
fp=open(sys.argv[1],'r')
for line in fp:
    arr=line.strip().split('\t')
    tnum+=int(arr[0])
    all[arr[1]]=int(arr[0])
fp.close()


snum=0
ind={}
fp=open(sys.argv[2],'r')
for line in fp:
    arr=line.strip().split('\t')
    snum+=int(arr[0])
    ind[arr[1]]=int(arr[0])
fp.close()

print 'Fisher-test:'
for id in ind:
    table=[[],[]]
    table[0].append(ind[id])
    table[0].append(all[id])
    table[1].append(snum-ind[id])
    table[1].append(tnum-all[id])
    oddsratio, pvalue = stats.fisher_exact(table)
    print str(ind[id])+'\t'+id+'\t'+str(oddsratio)+'\t'+str(pvalue)

print '\nBinormal-test:'
ratio=snum*1.0/tnum
for id in ind:
    pval = stats.binom.sf(ind[id],all[id],ratio)
    #print ratio
    #print str(ind[id])+'\t'+id+'\t'+str(pval)
    print id+'('+str(ind[id])+', '+str(pval)+')'
