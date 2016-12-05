
import subprocess


cmda = '''
time pydair parse -q data/simdata.%s.fa \\
                  -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \\
                  --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \\
                  -o results/simexe.%s \\
'''


cmdb = '''
time pydair stats -i results/simexe.%s.vdj.pydair -n simexestats%s\\
                  -o results/simexestats%s --contain_ambiguous_D \\
                  --all_seq --estimate-vdj-combination
'''


subset = [1000, 5000, 10000, 20000, 40000, 60000, 80000]
subset = [100000]
for s in subset:
    pydair_cmd = cmda % (s, s)
    pydair_cmd.rstrip()
    proc_sig = subprocess.call(pydair_cmd, shell = True)


for s in subset:
    pydair_cmd = cmdb % (s, s, s)
    pydair_cmd.rstrip()
    proc_sig = subprocess.call(pydair_cmd, shell = True)





