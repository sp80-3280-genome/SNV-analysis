"""
This script takes vcf and gff3 file as input; outputs SNP annotations.
By Hui Guo
5/16/2012
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


PROMOTER_LENGTH=2000
GENE='gene'
UTR5='five_prime_UTR'
UTR3='three_prime_UTR'
CDS='CDS'

chr2gene={} # list of (gid, strand, start, end)
chr2splice={} # list of position of splice sites
chr2promoter={} # list of coords (start, end)
chr2utr5={} # list of coords (start, end)
chr2utr3={} # list of coords (start, end)
chr2cds={} # list of (phase, start, end)

#print 'Reading gff3...'
fp=open('../../sc04_sc.mlc.cns.sgl.importdb.scga7.gff')
gff3=fp.readlines()
fp.close()
gkey=''
softwares=['EVM','PASA']
for line in gff3:
    if line.startswith('#') or line =='\n': continue
    arr=line.strip().split('\t')
    if arr[1] not in softwares: continue
    if arr[2] == GENE:
        gid=(arr[8].split(';'))[0][3:]
        gid=gid[:4]+'model'+gid[6:]
        gkey=arr[0]+'_'+gid
        if arr[0] not in chr2gene: chr2gene[arr[0]]=[(gid, arr[6], int(arr[3]), int(arr[4]))]
        else: chr2gene[arr[0]].append((gid, arr[6], int(arr[3]), int(arr[4])))
        if arr[6]=='+':
            pstart=int(arr[3])-PROMOTER_LENGTH
            pend=int(arr[3])-1
            chr2promoter[gkey]=(pstart,pend)
        else:
            pstart=int(arr[4])+1
            pend=int(arr[4])+PROMOTER_LENGTH
            chr2promoter[gkey]=(pstart,pend)

    elif arr[2] == UTR5:
        gid=(arr[8].split(';'))[0][3:]
        gkey=arr[0]+'_'+gid
        if gkey not in chr2utr5: chr2utr5[gkey]=[(int(arr[3]), int(arr[4]))]
        else: chr2utr5[gkey].append((int(arr[3]), int(arr[4])))
        if arr[6] == '+':
            if gkey not in chr2splice: chr2splice[gkey]=[int(arr[4])+1, int(arr[4])+2]
            else: 
                chr2splice[gkey].append(int(arr[4])+1)
                chr2splice[gkey].append(int(arr[4])+2)
        else:
            if gkey not in chr2splice: chr2splice[gkey]=[int(arr[3])-1, int(arr[3])-2]
            else: 
                chr2splice[gkey].append(int(arr[3])-1)
                chr2splice[gkey].append(int(arr[3])-2)
    elif arr[2] == UTR3:
        gid=(arr[8].split(';'))[0][3:]
        gkey=arr[0]+'_'+gid
        if gkey not in chr2utr3: chr2utr3[gkey]=[(int(arr[3]), int(arr[4]))]
        else: chr2utr3[gkey].append((int(arr[3]), int(arr[4])))
        if arr[6] == '+':
            if gkey not in chr2splice: chr2splice[gkey]=[int(arr[3])-1, int(arr[3])-2]
            else:
                chr2splice[gkey].append(int(arr[3])-1)
                chr2splice[gkey].append(int(arr[3])-2)
        else:
            if gkey not in chr2splice: chr2splice[gkey]=[int(arr[4])+1, int(arr[4])+2]
            else: 
                chr2splice[gkey].append(int(arr[4])+1)
                chr2splice[gkey].append(int(arr[4])+2)
    elif arr[2] == CDS:
        gid=(arr[8].split(';'))[0][3:]
        gkey=arr[0]+'_'+gid[4:]
        if gkey not in chr2cds: chr2cds[gkey]=[(arr[7], int(arr[3]), int(arr[4]))]
        else: chr2cds[gkey].append((arr[7], int(arr[3]), int(arr[4])))
        if gkey not in chr2splice: chr2splice[gkey]=[int(arr[3])-1, int(arr[3])-2, int(arr[4])+1, int(arr[4])+2]
        else:
            chr2splice[gkey].append(int(arr[3])-1)
            chr2splice[gkey].append(int(arr[3])-2)
            chr2splice[gkey].append(int(arr[4])+1)
            chr2splice[gkey].append(int(arr[4])+2)

#print 'Reading Genome sequences...'
chr2seq={}
fp=open('../../sc03_sc.mlc.cns.sgl.utg.scga7.importdb.fa')
for seq_record in SeqIO.parse(fp, "fasta"):
    chr2seq[seq_record.id] = str(seq_record.seq).upper()
fp.close()


def Warn():
    print 'ThisRefn not match the nucleiotide on genome sequence... Skipped...'

def InGenes(Chr, pos):
    if Chr not in chr2gene: return 0,0,0,0
    for g, d, s, e in chr2gene[Chr]:
        if s<= pos <=e:
            return g, d, s, e
    return 0, 0, 0, 0

#print 'Parsing vcf and writing results...'
fp=open(sys.argv[1])
vcf=fp.readlines()
fp.close()
INDEL=0
for line in vcf:
    if line.startswith('#'): continue
    arr=line.split('\t')
    if arr[7].startswith('INDEL'): INDEL=1
    ThisChr=arr[0]
    ThisPos=int(arr[1])
    ThisRefn=arr[3]
    ThisAltn=arr[4]
    ThisGid='.'
    ThisStrand='+'
    ThisType='Inter-genic'
    ThisPhase='.'
    ThisRefaa='.'
    ThisAltaa='.'
    ThisGenePos='.'
    ThisCDSPos='.'
    ThisProteinPos='.'
    
    if INDEL:
        reflen=len(ThisRefn)
        altlen=len(ThisAltn)
        if reflen > altlen:
            ThisAltn='-'+ThisRefn[altlen:]
            ThisPos=ThisPos+altlen-1
        else:
            ThisAltn='+'+ThisAltn[reflen:]
            ThisPos=ThisPos+reflen-1
    
    gid, strand, start, end = InGenes(ThisChr, ThisPos)
    if gid:
        key=ThisChr+'_'+gid
        ThisStrand=strand
        ThisGid=gid
        
        if key in chr2splice:
            if ThisPos in chr2splice[key]:
                ThisType='Splice'
                if ThisStrand == '+': ThisGenePos = ThisPos-start+1
                else: ThisGenePos = end-ThisPos+1
        
        if key in chr2promoter:    
            s, e = chr2promoter[key]
            if s<= ThisPos <=e:
                ThisType='Promoter'
                if ThisStrand == '+': ThisGenePos = -(e-ThisPos+1)
                else: ThisGenePos = -(ThisPos-s+1)
        
        if key in chr2utr5:
            chr2utr5[key]=sorted(chr2utr5[key], key=lambda x: x[0])
            for s,e in chr2utr5[key]:
                if s<= ThisPos <=e:
                    ThisType="5-UTR"
                    ThisCDSPos = e-s+1
                  
                    if ThisStrand == '+': 
                        ThisGenePos = ThisPos-start+1
                        ThisPhase = ThisPos-s+1 
                    else: 
                        ThisGenePos = end-ThisPos+1
                        ThisPhase = e-ThisPos+1
                    break
        
        if key in chr2utr3:
            chr2utr3[key]=sorted(chr2utr3[key], key=lambda x: x[0])
            for s,e in chr2utr3[key]:
                if s<= ThisPos <=e:
                    ThisType="3-UTR"
                    ThisCDSPos = e-s+1
                
                    if ThisStrand == '+': 
                        ThisGenePos = ThisPos-start+1
                        ThisPhase = ThisPos-s+1
                    else: 
                        ThisGenePos = end-ThisPos+1
                        ThisPhase = e-ThisPos+1
                    break

        count=0
        if key not in chr2cds: continue
        chr2cds[key]=sorted(chr2cds[key], key=lambda x: x[1])
        
        #print chr2cds['Chr01_Goraiv21000007m']
        #mycds=''
        #mypep=''
        #for myp,mys,mye in chr2cds['Chr01_Goraiv21000007m']:
        #    myexon=chr2seq['Chr01'][mys-1:mye]
        #    myexon=str(Seq(myexon, IUPAC.unambiguous_dna).reverse_complement())
        #    mycds+=myexon
        #mypep=str(Seq(mycds, IUPAC.unambiguous_dna).translate())
        #print mycds
        #print mypep
        #sys.exit()

        for p,s,e in chr2cds[key]:
            count+=1
            p=int(p)
            if s<= ThisPos <=e:
                ThisType="CDS"
                
                frontExonLen=0
                altn=ThisAltn.split(',')
                altcodon=[]
                altaa=[]
                if ThisStrand == '+': 
                    ThisGenePos = ThisPos-start+1
                    
                    m = (ThisPos-s-p+1)%3
                    if m==0: ThisPhase=1
                    elif m==1: ThisPhase=0
                    elif m==2: ThisPhase=2
                    
                    if INDEL==0:
                        if ThisPhase == 0:
                            Refcodon=chr2seq[ThisChr][ThisPos-1:ThisPos+2]
                            for a in altn:
                                altcodon.append(a+chr2seq[ThisChr][ThisPos]+chr2seq[ThisChr][ThisPos+1])
                        elif ThisPhase == 1:
                            Refcodon=chr2seq[ThisChr][ThisPos-3:ThisPos]
                            for a in altn:
                                altcodon.append(chr2seq[ThisChr][ThisPos-3]+chr2seq[ThisChr][ThisPos-2]+a)
                        elif ThisPhase == 2:
                            Refcodon=chr2seq[ThisChr][ThisPos-2:ThisPos+1]
                            for a in altn:
                                altcodon.append(chr2seq[ThisChr][ThisPos-2]+a+chr2seq[ThisChr][ThisPos])
                        if chr2seq[ThisChr][ThisPos-1] != ThisRefn: Warn()
                        else:
                            try:
                                ThisRefaa=str(Seq(Refcodon, IUPAC.unambiguous_dna).translate())
                            except:
                                ThisRefaa="-"
                            for codon in altcodon:
                                try:
                                    altaa.append(str(Seq(codon, IUPAC.unambiguous_dna).translate()))
                                except:
                                    altaa.append("-")
                    
                    for i in range(count-1):
                        pp,ss,ee=chr2cds[key][i]
                        thisExonLen=ee-ss+1
                        frontExonLen+=thisExonLen
                    ThisCDSPos=frontExonLen+ThisPos-s+1

                else: 
                    ThisGenePos = abs(ThisPos-end)+1
                    
                    m = (e-ThisPos-p+1)%3
                    if m==2: ThisPhase=2
                    elif m==1: ThisPhase=0
                    elif m==0: ThisPhase=1
                    
                    if INDEL==0:
                        if ThisPhase == 0:
                            Refcodon=chr2seq[ThisChr][ThisPos-3:ThisPos]
                            for a in altn:
                                altcodon.append(chr2seq[ThisChr][ThisPos-3]+chr2seq[ThisChr][ThisPos-2]+a)
                        elif ThisPhase == 1:
                            Refcodon=chr2seq[ThisChr][ThisPos-1:ThisPos+2]
                            for a in altn:
                                altcodon.append(a+chr2seq[ThisChr][ThisPos]+chr2seq[ThisChr][ThisPos+1])
                        elif ThisPhase == 2:
                            Refcodon=chr2seq[ThisChr][ThisPos-2:ThisPos+1]
                            for a in altn:
                                altcodon.append(chr2seq[ThisChr][ThisPos-2]+a+chr2seq[ThisChr][ThisPos])
                    
                        Refcodon=str(Seq(Refcodon, IUPAC.unambiguous_dna).reverse_complement())
                        altcodonRev=[]
                        for c in altcodon:
                            altcodonRev.append(str(Seq(c, IUPAC.unambiguous_dna).reverse_complement()))
                        if chr2seq[ThisChr][ThisPos-1] != ThisRefn: Warn()
                        else:
                            try:
                                ThisRefaa=str(Seq(Refcodon, IUPAC.unambiguous_dna).translate())
                            except:
                                ThisRefaa="-"
                            for codon in altcodonRev:
                                try:
                                    altaa.append(str(Seq(codon, IUPAC.unambiguous_dna).translate()))
                                except:
                                    altaa.append("-")
                    
                    for i in range(count,len(chr2cds[key])):
                        pp,ss,ee=chr2cds[key][i]
                        thisExonLen=ee-ss+1
                        frontExonLen+=thisExonLen
                    ThisCDSPos=frontExonLen+e-ThisPos+1
                        
                if INDEL==0: ThisAltaa=','.join(altaa)
                else: ThisAltaa='Indel'
                
                proteinPos=ThisCDSPos%3
                if proteinPos == 0: ThisProteinPos=ThisCDSPos/3
                else: ThisProteinPos=ThisCDSPos/3+1
                break

        # Intron
        if ThisGenePos =='.':
            tmparr=[(y,z) for x,y,z in chr2cds[key]]
            if key not in chr2utr3 and key not in chr2utr5:
                Exons=tmparr
            elif key not in chr2utr3:
                Exons=chr2utr5[key]+tmparr
            elif key not in chr2utr5:
                Exons=tmparr+chr2utr3[key]
            else:
                Exons=chr2utr5[key]+tmparr+chr2utr3[key]
            if ThisStrand == '+': 
                ThisGenePos = ThisPos-start+1
                for j in range(len(Exons)-1):
                    s1, e1 = Exons[j]
                    s2, e2 = Exons[j+1]
                    if e1< ThisPos <s2:
                        ThisCDSPos=s2-e1+1 # length of intron
                        ThisPhase=ThisPos-e1+1
            else: 
                ThisGenePos = abs(ThisPos-end)+1
                for j in range(len(Exons)-1):
                    s1, e1 = Exons[j]
                    s2, e2 = Exons[j+1]
                    if e2< ThisPos <s1:
                        ThisCDSPos=s1-e2+1 # length of intron
                        ThisPhase=s1-ThisPos+1
            ThisType="Intron"
    if INDEL!=0: ThisAltaa='Indel'     
    INDEL=0
    ThisLine=[ThisChr, str(ThisPos), ThisStrand, ThisRefn, ThisAltn,\
              ThisGid, ThisType, str(ThisPhase), str(ThisGenePos),\
              str(ThisCDSPos), str(ThisProteinPos), ThisRefaa, ThisAltaa]
    print '\t'.join(ThisLine)
