
import subprocess


cmd = '''
pydair parse -q data/simdata.sub.fa \\
             -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \\
             --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \\
             --v-wordsize %s --v-match-score %s --v-mismatch-score %s --v-gap-open-penalty %s --v-gap-extend-penalty %s \\
             --j-wordsize %s --j-match-score %s --j-mismatch-score %s --j-gap-open-penalty %s --j-gap-extend-penalty %s \\
             -o results/sim%s

'''

v = [{'wordsize': 10, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 6},
     {'wordsize': 10, 'match': 3, 'mismatch': -3, 'gapopen': 3, 'gapextend': 6},
     {'wordsize': 10, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 3},
     {'wordsize': 10, 'match': 2, 'mismatch': -7, 'gapopen': 6, 'gapextend': 6},
     {'wordsize': 20, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 6},
     {'wordsize': 30, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 6}]
j = [{'wordsize':  7, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 6},
     {'wordsize':  7, 'match': 3, 'mismatch': -3, 'gapopen': 3, 'gapextend': 6},
     {'wordsize':  7, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 3},
     {'wordsize':  7, 'match': 2, 'mismatch': -7, 'gapopen': 6, 'gapextend': 6},
     {'wordsize':  4, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 6},
     {'wordsize': 20, 'match': 3, 'mismatch': -3, 'gapopen': 6, 'gapextend': 6}]


for vi in range(len(v)):
    for ji in range(len(j)):
        if ji != 4:
            continue
        prefix = str(vi) + '_' + str(ji)
        pydair_cmd = cmd % (v[vi]['wordsize'], v[vi]['match'], v[vi]['mismatch'], v[vi]['gapopen'], v[vi]['gapextend'],
                            j[ji]['wordsize'], j[ji]['match'], j[ji]['mismatch'], j[ji]['gapopen'], j[ji]['gapextend'],
                            prefix)
        pydair_cmd.rstrip()
        print pydair_cmd
        proc_sig = subprocess.call(pydair_cmd, shell = True)
        
        
        


