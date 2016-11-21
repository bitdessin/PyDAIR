import sys
import os
import re

input_file = sys.argv[1]
output_file = sys.argv[2]


with open(output_file, 'w') as fho:
    with open(input_file, 'r') as fh:
        for buf in fh:
            fields = buf.strip().split('\t')
            if 'Seq-' not in fields[0]:
                continue
            seq_id = fields[0]
            seq    = fields[1]
            vgene  = fields[21]
            dgene  = fields[22]
            jgene  = fields[23]
            cdr3   = fields[71].replace(' ', '')
            warns  = fields[74]
            
            if 'out of frame' in warns:
                stop_codon = 'has_stop_codon'
            else:
                stop_codon = 'without_codon'
            fa_data = '>' + seq_id + '; ' + vgene + '; ' + dgene + '; ' + jgene + '; ' + cdr3 + '; ' + stop_codon + '\n'
            fa_data = fa_data + seq + '\n'
            fho.write(fa_data)

