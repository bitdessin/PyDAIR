#!/usr/bin/env python
import os
import sys
import argparse

log_file = sys.argv[1]
in_fasta = sys.argv[2]
out_fasta = sys.argv[3]


def read_cutadaptlog(logf, pdict):
    with open(logf, 'r') as infh:
        for buf in infh:
            rd = buf.split('\t')
            if len(rd) > 7:
                rd2 = rd[0].split(' ')
                seqid = rd2[0]
                primer = rd[7].strip()
                pdict[seqid] = primer
    return pdict
    

def printout_fq(inf, outf, pdict):
    igm = outf + '.igm.fq'
    igz = outf + '.igz.fq'
    
    igmfh = open(igm, 'w')
    igzfh = open(igz, 'w')
    
    with open(inf, 'r') as fqin:
        cnt = 0
        seqid  = None
        seqdat = None
        for buf in fqin:
            cnt += 1
            if cnt % 4 == 1:
                if seqdat is not None:
                    if seqid in pdict:
                        if pdict[seqid] == 'IGM':
                            igmfh.write(seqdat)
                        if pdict[seqid] == 'IGZ':
                            igzfh.write(seqdat)
                seqid  = buf.split(' ')[0][1:]
                seqdat = buf
            else:
                seqdat = seqdat + buf
    
    igmfh.close()
    igzfh.close()


    
pdict = {}
pdict = read_cutadaptlog(log_file, pdict)
printout_fq(in_fasta, out_fasta, pdict)

    
    
