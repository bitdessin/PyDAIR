import os
import sys
import argparse
import subprocess
import unittest
from PyDAIR.seq.IgSeq import IgSeq
from PyDAIR.io.PyDAIRIO import *
from PyDAIR.utils.PyDAIRUtils import *
from PyDAIR.utils.PyDAIRArgs import *
from PyDAIR.app.PyDAIRAPP import *

_data_path = os.path.join(os.path.dirname(__file__), 'data/samples')
_db_path = os.path.join(os.path.dirname(__file__), 'data/db')
_result_path = os.path.join(os.path.dirname(__file__), 'data/results')


class Test_bin(unittest.TestCase):
    
    def setUp(self):
        pass
    
    
    def test_bin_parseseq(self):
        cmd = 'pydair parse -s fugu -q ' + _data_path + '/sample.1.fa -o ' + _result_path + '/test_output_bin_parseseq_ -f pydair '
        cmd += '-v ' + _db_path + '/v.fa -d ' + _db_path + '/d.fa -j ' + _db_path + '/j.fa '
        cmd += '--v-blastdb ' + _db_path + '/v '
        cmd += '--v-match-score 3 --v-mismatch-score -3 '
        cmd += '--v-gap-open-penalty 6 --v-gap-extend-penalty 6 '
        cmd += '--v-wordsize 21 --v-evalue-cutoff 1e-60 '
        cmd += '--d-blastdb ' + _db_path + '/d '
        cmd += '--d-match-score 1 --d-mismatch-score -1 '
        cmd += '--d-gap-open-penalty 0 --d-gap-extend-penalty 2 '
        cmd += '--d-wordsize 4 --d-evalue-cutoff 1 '
        cmd += '--j-blastdb ' + _db_path + '/j '
        cmd += '--j-match-score 3 --j-mismatch-score -3 '
        cmd += '--j-gap-open-penalty 6 --j-gap-extend-penalty 6 '
        cmd += '--j-wordsize 7 --j-evalue-cutoff 1e-5 '
        print(cmd)
        subprocess.call(cmd, shell = True)
    
    
    
    def test_bin_stats(self):
        cmd = 'pydair stats -i ' + _data_path + '/sample.1.pydair ' + _data_path + '/sample.2.pydair ' + _data_path + '/sample.3.pydair '
        cmd += '-n sample_1 smaple_2 sample_3 '
        cmd += '-o ' + _result_path + '/test_output_bin_analysis '
        print(cmd)
        subprocess.call(cmd, shell = True)
    
    
    
    def test_bin_stats_2(self):
        cmd = 'pydair stats -i ' + _data_path + '/sample.1.pydair ' + _data_path + '/sample.2.pydair ' + _data_path + '/sample.3.pydair '
        cmd += '-n sample_1 smaple_2 sample_3 '
        cmd += '-o ' + _result_path + '/test_output_bin_analysis_hasambigoD '
        cmd += '--contain_ambiguous_D '
        cmd += '--figure-style ggplot '
        print(cmd)
        subprocess.call(cmd, shell = True)





if __name__ == '__main__':
    unittest.main()

