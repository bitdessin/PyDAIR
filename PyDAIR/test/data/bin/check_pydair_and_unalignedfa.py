"""
This Python script is used to check the sequences which are
consisted in PYDAIR format file but are not consisted in
un-aligned FASTA file.

Sequences satisfied the following two conditions simultaneously
are consisted in PYDAIR format but not in un-aligned FASTA file.

- the V-aligned region and J-aligned region have overlap
- both YYC and WGxG are cannot be found


python test_bin.py
python ./data/bin/check_pydair_and_unalignedfa.py
"""



unaligned_f = './data/results/test_output_bin_parseseq_.unaligned.fa'
pydair_f = './data/results/test_output_bin_parseseq_.vdj.pydair'

seqid = {}

with open(unaligned_f, 'r') as fafh:
    for buf in fafh:
        if buf[0:1] == '>':
            seqid[buf[1:].split(' ')[0]] = 1


with open(pydair_f, 'r') as pyfh:
    for buf in pyfh:
        if buf[0:2] == 'QN':
            sid = buf[3:].rsplit()[0]
            if sid not in seqid:
                print sid




