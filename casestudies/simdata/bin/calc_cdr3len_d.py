import sys
import os
import re

sim_file = sys.argv[1]
res_file = sys.argv[2]
fa_file  = sys.argv[3]
out_file = sys.argv[4]

fa_dict = {}


def check(true = None, estimated = None):
    tag = ''
    if estimated == '.':
        tag = 'ambiguous'
    elif true in estimated or estimated in true:
        tag = 'correct'
    else:
        tag = 'incorrect'
    return tag



with open(res_file, 'r') as fh:
    for buf in fh:
        fields  = buf.strip().split('\t')
        seq_id  = fields[0].replace(';', '')
        dgene   = fields[4]
        fa_dict[seq_id] = {'d': dgene, 'cdr3len': -1}


with open(fa_file, 'r') as fh:
    seq_id = None
    for buf in fh:
        if buf[0:1] == '>':
            fields = buf.strip().split(' ')
            seq_id = fields[0].replace(';', '').replace('>', '')
        else:
            buf = buf.strip()
            if seq_id in fa_dict:
                fa_dict[seq_id]['cdr3len'] = len(buf)
            else:
                fa_dict[seq_id] = {'d': None, 'cdr3len': len(buf)}


outfh = open(out_file, 'w')

with open(sim_file, 'r') as fh:
    for buf in fh:
        fields = buf.strip().split('\t')
        if 'Seq-' not in fields[0]:
            continue
        seq_id = fields[0]
        dgene  = fields[22]
        
        # check estimation of stop codon stats.
        tag = 'ambiguous'
        if seq_id in fa_dict:
            tag = check(dgene, fa_dict[seq_id]['d'])
            cdr3len = fa_dict[seq_id]['cdr3len']
        
        outfh.write(seq_id + '\t' + tag + '\t' + str(cdr3len) + '\n')
        

'''
library(ggplot2)

x <- read.table('joinsim.sub.dstats.txt')
x[, 2] <- factor(x[, 2], levels = c('correct', 'incorrect', 'ambiguous'))
x <- x[x[, 3] < 100, ]
colnames(x) <- c('id', 'Result', 'Length')
g <- ggplot(x, aes(x = Length, fill = Result))
g <- g + geom_histogram(position = 'fill', binwidth = 3)
g <- g + xlab('Sequence length of unaligned region')
g <- g + ylab('Probability')
g <- g + theme_bw() + scale_fill_manual(values = c('#FC8D59', '#FFFFBF', '#91CF60'))
pdf('figure4.pdf', width = 6, height = 4)
plot(g)
dev.off()

'''





