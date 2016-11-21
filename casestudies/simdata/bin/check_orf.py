import sys
import os
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

sim_file = sys.argv[1]
out_file = sys.argv[2]


re_start_codon = re.compile('M')
re_stop_codons = '*'


with open(sim_file, 'r') as fh:
    fields = buf.strip().split('\t')
        if 'Seq-' not in fields[0]:
            continue
    seq_id = fields[0]
    seq_nl = fields[1]
    
    if seq_id in target_ids:
        orf_i = None
        
        for orf in range(3):
            orf_i = None
            seq_nl_frame = seq_nl[(orf):int(math.floor((len(seq_nl) - orf) / 3) * 3 + orf)]
            seq_aa       = str(Seq(seq_nl_frame, generic_dna).translate())
            
            for m in re_start_codon.finditer(seq_aa):
                seq_aa_orf = seq_aa[(m.start()):]
                if re_stop_codons not in seq_aa_orf:
                    orf_i = i
                    break
                else:
                    continue
            
            






