import os
import math
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from PyDAIR.seq.IgSeq import IgSeq
from PyDAIR.io.PyDAIRIO import *
from PyDAIR.stats.PyDAIRStats import *

_data_path = os.path.join(os.path.dirname(__file__), 'data/samples')
_result_path = os.path.join(os.path.dirname(__file__), 'data/results')

class Test_pydair_stats(unittest.TestCase):
    
    
    def setUp(self):
        self.pydair_input_files = [_data_path + '/sample.1.pydair',
                                   _data_path + '/sample.2.pydair',
                                   _data_path + '/sample.3.pydair']
        self.pydair_id = ['Sample 1', 'Sample 2', 'Sample 3']
        self.pydair_div_f = [_result_path + '/test_output_stat_f.1.txt',
                             _result_path + '/test_output_stat_f.2.txt',
                             _result_path + '/test_output_stat_f.3.txt']
        self.pydair_div_s = [_result_path + '/test_output_stat_s.1.txt',
                             _result_path + '/test_output_stat_s.2.txt',
                             _result_path + '/test_output_stat_s.3.txt']
    
    
    
    
    def test_stats_freq(self):
        print('test_stats_freq')
        pydair_files = self.pydair_input_files
        pydair_id    = self.pydair_id
        bstats = PyDAIRStats(pydair_files, pydair_id)
        
        # here should be None object
        for bsample in bstats.samples:
            print(bsample.get_summary('vdj_rarefaction'))
        print(bstats.get_summary('vdj_rarefaction'))
        
        bstats.rarefaction_study('vdj', 2)
        
        print(len(bstats.samples))
        print(bstats.samples.len())
        for bsample in bstats.samples:
            print(bsample.vdj.head(5))
            print(bsample.cdr3.head(5))
            print(bsample.indels.head(5))
            
            print(bsample.get_summary('v').head(5))
            print(bsample.get_summary('d').head(5))
            print(bsample.get_summary('j').head(5))
            print(bsample.get_summary('vdj').head(5))
            print(bsample.get_summary('cdr3_nucl_len').head(5))
            print(bsample.get_summary('cdr3_prot_len').head(5))
            print(bsample.get_summary('v_del_len').head(5))
            print(bsample.get_summary('j_del_len').head(5))
            print(bsample.get_summary('vj_ins_len').head(5))
            print(bsample.get_summary('vdj_rarefaction').head(5))
        
        print(bstats.get_summary('v').head(5))
        print(bstats.get_summary('d').head(5))
        print(bstats.get_summary('j').head(5))
        print(bstats.get_summary('vdj').head(5))
        print(bstats.get_summary('cdr3_nucl_len').head(5))
        print(bstats.get_summary('cdr3_prot_len').head(5))
        print(bstats.get_summary('v_del_len').head(5))
        print(bstats.get_summary('j_del_len').head(5))
        print(bstats.get_summary('vj_ins_len').head(5))
        print(bstats.get_summary('vdj_rarefaction').head(5))
        
        

if __name__ == '__main__':
    unittest.main()



