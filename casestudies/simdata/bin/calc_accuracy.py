import sys
import os
import re

sim_file = sys.argv[1]
res_file = sys.argv[2]
out_file = sys.argv[3]

fa_dict = {}

vstat    = {'correct': 0, 'mistake': 0, 'fault': 0, 'ambigo': 0}
dstat    = {'correct': 0, 'mistake': 0, 'fault': 0, 'ambigo': 0}
jstat    = {'correct': 0, 'mistake': 0, 'fault': 0, 'ambigo': 0}
cdr3stat = {'correct': 0, 'mistake': 0, 'fault': 0, 'ambigo': 0}
scstat   = {'TP': 0, 'TN': 0, 'FP': 0, 'FN': 0, 'fault': 0}



def check(d, true = None, estimated = None, t = 'seq', seqid = None):
    if t == 'seq':
        if estimated == '.':
            d['ambigo'] += 1
        elif true in estimated or estimated in true:
            d['correct'] += 1
            print 'correct: ' + seq_id
        else:
            d['mistake'] += 1
            print 'mistake: ' + seq_id
            print 'true:    ' + true
            print 'esti:    ' + estimated
    
    elif t == 'orf':
        if true == 'N' and estimated == 'N':
            d['TP'] += 1
            print 'TP: ' + seq_id
        elif true == 'N' and estimated == '*':
            d['FN'] += 1
            print 'FN: ' + seq_id
        elif true == '*' and estimated == 'N':
            d['FP'] += 1
            print 'FP: ' + seq_id
        elif true == '*' and estimated == '*':
            d['TN'] += 1
            print 'TN: ' + seq_id
    
    return d
    
    

with open(res_file, 'r') as fh:
    for buf in fh:
        fields  = buf.strip().split('\t')
        seq_id  = fields[0].replace(';', '')
        orf     = fields[1]
        orfcode = fields[2]
        vgene   = fields[3]
        dgene   = fields[4]
        jgene   = fields[5]
        cdr3    = fields[7]
        if orf == '.':
            orf = '*'
        else:
            orf = 'N'
        fa_dict[seq_id] = {'v': vgene, 'd': dgene, 'j': jgene, 'cdr3': cdr3, 'stop': orf}


with open(sim_file, 'r') as fh:
    for buf in fh:
        fields = buf.strip().split('\t')
        if 'Seq-' not in fields[0]:
            continue
        seq_id = fields[0]
        vgene  = fields[21]
        dgene  = fields[22]
        jgene  = fields[23]
        cdr3   = fields[71].replace(' ', '')
        stop   = fields[74]
        if stop == '[]':
            stop = 'N'
        else:
            stop = '*'
        
        # check estimation of stop codon stats.
        if seq_id in fa_dict:
            vstat    = check(vstat, vgene, fa_dict[seq_id]['v'], 'seq')
            dstat    = check(dstat, dgene, fa_dict[seq_id]['d'], 'seq')
            jstat    = check(jstat, jgene, fa_dict[seq_id]['j'], 'seq')
            cdr3stat = check(cdr3stat, cdr3, fa_dict[seq_id]['cdr3'], 'seq')
            scstat   = check(scstat, stop, fa_dict[seq_id]['stop'], 'orf', seq_id)
        else:
            print "Fault: " + seq_id
            vstat['fault']    += 1
            dstat['fault']    += 1
            jstat['fault']    += 1
            cdr3stat['fault'] += 1
            scstat['fault']   += 1


with open(out_file, 'w') as fh:
    fh.write('   \t  TP  \t  TN  \t  FP  \t  FN  \t  fault\n')
    fh.write('ORF\t' + str(scstat['TP']) + '\t' + str(scstat['TN']) + '\t' + \
             str(scstat['FP']) + '\t' + str(scstat['FN']) + '\t' + str(scstat['fault']) + '\n\n')
    fh.write('feature\tcorrect\tmistake\tambigo\tfault\n')
    fh.write('V      \t' + str(vstat['correct']) +'\t' + str(vstat['mistake']) + '\t' + str(vstat['ambigo']) + '\t' + str(vstat['fault']) + '\n')
    fh.write('D      \t' + str(dstat['correct']) +'\t' + str(dstat['mistake']) + '\t' + str(dstat['ambigo']) + '\t' + str(dstat['fault']) + '\n')
    fh.write('J      \t' + str(jstat['correct']) +'\t' + str(jstat['mistake']) + '\t' + str(jstat['ambigo']) + '\t' + str(jstat['fault']) + '\n')
    fh.write('CDR3   \t' + str(cdr3stat['correct']) +'\t' + str(cdr3stat['mistake']) + '\t' + str(cdr3stat['ambigo']) + '\t' + str(cdr3stat['fault']) + '\n')







