import sys
import os
import re

## python bin/calc_accuracy_details.py ./data/simdata.txt ./results/simdata.vdj.pydair.simple ./results/estperformance


sim_file = sys.argv[1]
res_file = sys.argv[2]
out_file = sys.argv[3]


v = set(['unidentifiable'])
d = set(['unidentifiable'])
j = set(['unidentifiable'])

vstat = {'truth': {}, 'sen': {}, 'fd': {}, 'details': {}}
dstat = {'truth': {}, 'sen': {}, 'fd': {}, 'details': {}}
jstat = {'truth': {}, 'sen': {}, 'fd': {}, 'details': {}}


est = {}
with open(res_file, 'r') as fh:
    for buf in fh:
        buf = buf.strip()
        if buf[0:6] == '#BEGIN':
            seq_id = ''
            orf   = ''
            vgene = ''
            dgene = ''
            jgene = ''
            cdr3  = ''
        if buf[0:2] == 'QN':
            seq_id = buf[3:].replace(';', '')
        if buf[0:2] == 'VN':
            vgene = buf[3:]
        if buf[0:2] == 'DN':
            dgene = buf[3:]
        if buf[0:2] == 'JN':
            jgene = buf[3:]
        if buf[0:2] == 'OP':
            orf = buf[3:]
            if orf == '.':
                orf = '*'
            else:
                orf = 'N'
        if buf[0:4] == '#END':
            if '|' in vgene:
                vgene = vgene.split('|')[1]
            if '|' in dgene:
                dgene = dgene.split('|')[1]
            if '|' in jgene:
                jgene = jgene.split('|')[1]
            est[seq_id] = {'v': vgene, 'd': dgene, 'j': jgene}
            v.update([vgene])
            d.update([dgene])
            j.update([jgene])



truth = {}
with open(sim_file, 'r') as fh:
    for buf in fh:
        fields  = buf.strip().split('\t')
        if 'Seq-' not in fields[0]:
            continue
        seq_id = fields[0].replace(';', '')
        vgene  = fields[21]
        dgene  = fields[22]
        jgene  = fields[23]
        
        truth[seq_id] = {'v': vgene, 'd': dgene, 'j': jgene}
        v.update([vgene])
        d.update([dgene])
        j.update([jgene])



def write_stats(gvec, g, truth, est, f):

    senvec = {}
    for seq_id in truth.keys():
        truth_g = truth[seq_id][g]
        if seq_id in est:
            est_g = est[seq_id][g]
        else:
            est_g = 'unidentifiable'
        
        if truth_g not in senvec:
            senvec[truth_g] = {'truth': 0, 'correctly': 0}
        senvec[truth_g]['truth'] += 1
        if truth_g == est_g:
            senvec[truth_g]['correctly'] += 1
    
    fdvec = {}
    for seq_id in est.keys():
        truth_g = truth[seq_id][g]
        est_g   = est[seq_id][g]
        
        if est_g not in fdvec:
            fdvec[est_g] = {'estimated': 0, 'correctly': 0, 'incorrectly': 0, 'details': {}}
        fdvec[est_g]['estimated'] += 1
        if est_g == truth_g:
            fdvec[est_g]['correctly'] += 1
        else:
            fdvec[est_g]['incorrectly'] += 1
            if truth_g not in fdvec[est_g]['details']:
                fdvec[est_g]['details'][truth_g] = 0
            fdvec[est_g]['details'][truth_g] += 1

    with open(f, 'w') as fh:
        for vname in sorted(gvec):
            txt = ''
            if vname in fdvec:
                for kk, vv in fdvec[vname]['details'].iteritems():
                    txt += kk + ':' + str(vv) + ';'
            
            if vname in senvec:
                txt2 = vname + '\t' + str(senvec[vname]['truth'])  + '\t' + str(senvec[vname]['correctly']) + '\t'
            else:
                txt2 = vname + '\t\t\t'
            if vname in fdvec:
                txt2 = txt2 + str(fdvec[vname]['estimated']) + '\t' + str(fdvec[vname]['correctly'])
                txt2 = txt2 + '\t' + str(fdvec[vname]['incorrectly']) + '\t' + txt + '\n'
            else:
                txt2 = txt2 + '\n'
            fh.write(txt2)
            

write_stats(v, 'v', truth, est, out_file + '.v.csv')
write_stats(d, 'd', truth, est, out_file + '.d.csv')
write_stats(j, 'j', truth, est, out_file + '.j.csv')
    
    







