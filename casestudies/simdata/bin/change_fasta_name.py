import os
import sys

f = sys.argv[1]
o = f + '.new'

with open(o, 'w') as outfh:
    with open(f, 'r') as fh:
        for buf in fh:
            if buf[0:1] == '>':
                rd = buf.split('|')
                outfh.write('>' + rd[1] + '\n')
            else:
                outfh.write(buf)
    


