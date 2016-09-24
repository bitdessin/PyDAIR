import os
import math
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from PyDAIR.seq.IgSeq import IgSeq
from PyDAIR.io.PyDAIRIO import *
from PyDAIR.stats.PyDAIRStats import *

_data_path = os.path.join(os.path.dirname(__file__), 'data')
_result_path = os.path.join(os.path.dirname(__file__), 'data/results')

class Test_pydair_stats(unittest.TestCase):
    
    
    def setUp(self):
        # Outputs
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
    
    
    def test_stats_samplingresampling(self):
        pydair_files = self.pydair_input_files
        pydair_id    = self.pydair_id
        bstats = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        for bsample in bstats.samples:
            bsample.samplingresampling_study('vdj', 3)
            bsample.samplingresampling_study('cdr3', 3)
            print(bsample.div.samplingresampling['vdj'])
            print(bsample.div.samplingresampling['cdr3'])
    
    
    def test_stats_rarefaction(self):
        pydair_files = self.pydair_input_files
        pydair_id    = self.pydair_id
        bstats = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        for bsample in bstats.samples:
            bsample.rarefaction_study('vdj', 2)
            bsample.rarefaction_study('cdr3', 2)
            print(bsample.div.rarefaction['vdj'])
            print(bsample.div.rarefaction['cdr3'])
            
    
    def test_stats_freq(self):
        pydair_files = self.pydair_input_files
        pydair_id    = self.pydair_id
        bstats = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        print(len(bstats.samples))
        print(bstats.samples.len())
        for bsample in bstats.samples:
            print(bsample.vdj.head(5))
            print(bsample.get_freq('v'))
            print(bsample.get_freq('d'))
            print(bsample.get_freq('j'))
            print(bsample.get_freq('vdj').head(5))
    
    
    def test_cdr3_stats(self):
        pydair_files = self.pydair_input_files
        pydair_id    = self.pydair_id
        bstats = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        for bsample in bstats.samples:
            print(bsample.cdr3.head(5))
            print(bsample.get_freq('cdr3_nucl_len').head(5))
            print(bsample.get_freq('cdr3_prot_len').head(5))
    
    
    def test_stats_methods(self):
        pydair_files = self.pydair_input_files
        pydair_id    = self.pydair_id
        bstats = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        
        df_freq_v = bstats.get_freq('v')
        df_freq_d = bstats.get_freq('d')
        df_freq_j = bstats.get_freq('j', prob = True)
        df_freq_cdr3len_1 = bstats.get_cdr3len_freq(prob = True)
        df_freq_cdr3len_2 = bstats.get_cdr3len_freq(prob = True)
        
        print(df_freq_v)
        print(df_freq_d)
        print(df_freq_j)
        print(df_freq_cdr3len_1)
        print(df_freq_cdr3len_2)


if __name__ == '__main__':
    unittest.main()



