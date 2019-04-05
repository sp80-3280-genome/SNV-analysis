import sys


def uniq(x):
    iset=set(x)
    return list(iset)

C=[]
P=[]
F=[]

atgene2goterm={}
fp=open('ATH_GO_GOSLIM.txt', 'r')
for line in fp:
    arr=line.strip().split('\t')
    if arr[0] in atgene2goterm:
        if arr[7] in atgene2goterm[arr[0]]:
            atgene2goterm[arr[0]][arr[7]].append(arr[8])
        else:
            atgene2goterm[arr[0]][arr[7]]=[arr[8]]
    else:
        atgene2goterm[arr[0]]={}
        atgene2goterm[arr[0]][arr[7]]=[arr[8]]
    if arr[7]=='C' and arr[8] not in C: C.append(arr[8])
    if arr[7]=='P' and arr[8] not in P: P.append(arr[8])
    if arr[7]=='F' and arr[8] not in F: F.append(arr[8])
fp.close()

sb2at={}
#fp=open('Osativa_204_annotation_info.txt', 'r')
#fp=open('Sbicolor_79_annotation_info.txt', 'r')
fp=open('sc_at.blastn.filt.out', 'r')
for line in fp:
    arr=line.split('\t')
    sb2at[arr[0]]=arr[1][:-2]
fp.close()

sb2anno={}
for sbid in sb2at:
    atid=sb2at[sbid]
    if atid in atgene2goterm:
        sb2anno[sbid]=atgene2goterm[atid]

sbids=[]
fp=open(sys.argv[1], 'r')
for line in fp:
    arr=line.strip().split('\t')
    sbids.append(arr[0])
fp.close()

gotype = sys.argv[2]
if gotype == 'C': GO=C
elif gotype == 'F': GO=F
elif gotype == 'P': GO=P
else: sys.exit()

# rank annotation terms in each category
anno2num={}
for i in GO:
    for id in sbids:
        if id not in sb2anno: continue
        if gotype not in sb2anno[id]:continue
        if i in sb2anno[id][gotype]:
            #if i == 'other membranes': print id
            if i in anno2num:
                anno2num[i]+=1
            else:
                anno2num[i]=1
for i in anno2num:
    print str(anno2num[i])+'\t'+i
