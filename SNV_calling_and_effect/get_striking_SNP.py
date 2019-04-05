import sys, re

fp=open(sys.argv[1]) # read output from SNPAnnotationer
myfile=fp.readlines()
fp.close()

prefix=sys.argv[1][:-4]
frameShift=open(prefix+'_frameShift.txt','w')
initiationAlt=open(prefix+'_initiationAlt.txt','w')
spliceSite=open(prefix+'_spliceSite.txt','w')
stopCodon=open(prefix+'_stopCodon.txt','w')
stopExtension=open(prefix+'_stopExtension.txt','w')
for line in myfile:
    arr=line.strip().split('\t')
    if arr[6]=='CDS' and arr[-1]=='Indel':
        indellen=(len(arr[4])-1)%3
        if indellen != 0: frameShift.write(line)
    if arr[6]=='CDS' and arr[-3]=='1' and arr[-1]!=arr[-2] and arr[-2]=='M':
        initiationAlt.write(line)
    if arr[6]=='Splice' and re.search('[GTA]', arr[3]):
        spliceSite.write(line)
    if arr[-1]=='*' and arr[-2]!='*':
        stopCodon.write(line)
    if arr[-1]!='*' and arr[-2]=='*':
        stopExtension.write(line)
frameShift.close()
initiationAlt.close()
spliceSite.close()
stopCodon.close()
stopExtension.close()
