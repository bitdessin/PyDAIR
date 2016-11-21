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

class Test_app(unittest.TestCase):
    
    def setUp(self):
        pass
    
    
    def test_app_parseseq(self):
        # variables settings
        v_gene_align_args = PyDAIRBlastArgs(_db_path + '/v', 3, -3, 6, 6, 21, 1e-60)
        d_gene_align_args = PyDAIRBlastArgs(_db_path + '/d', 1, -1, 0, 2,  4, 1)
        j_gene_align_args = PyDAIRBlastArgs(_db_path + '/j', 3, -3, 6, 6,  7, 1e-5)
        q_fasta      = _data_path + '/sample.1.fa'
        v_gene_fasta = _db_path + '/v.fa'
        d_gene_fasta = _db_path + '/d.fa'
        j_gene_fasta = _db_path + '/j.fa'
        output_prefix = _result_path + '/test_output_app_parseseq'
        
        # PyDAIR arguemnts settings
        pydair_args = PyDAIRParseSeqArgs(q_fasta, v_gene_fasta, d_gene_fasta, j_gene_fasta,
                                       output_prefix, 'pydair',
                                       v_gene_align_args, d_gene_align_args, j_gene_align_args)
        
        # main processes
        pydairapp = PyDAIRAPPParseSeq(pydair_args)
        pydairapp.blast('v')
        pydairapp.blast('j')
        pydairapp.parse_VJ()
        pydairapp.write_pydair()
        pydairapp.write_fasta('unaligned_seq')
        pydairapp.blast('d')
        pydairapp.parse_VDJ()
        pydairapp.write_pydair()
        pydairapp.write_pydair(file_format = 'simple')
        pydairapp.write_fasta('cdr3', 'prot')
    
    
    
    def test_app_stats(self):
        sample_names = ['ID 1', 'ID 2', 'ID 3']
        pydair_files  = [_data_path + '/sample.1.pydair',
                         _data_path + '/sample.2.pydair',
                         _data_path + '/sample.3.pydair']
        pydair_args = PyDAIRStatsArgs(sample_names, pydair_files, True, False, False, 2,
                                      _result_path + '/test_output_app_analysis_hasambigoD', 'pdf')
        pydairapp = PyDAIRAPPStats(pydair_args)
        pydairapp.write_freq()
        pydairapp.write_cdr3_len_freq()
        pydair_args = PyDAIRStatsArgs(sample_names, pydair_files, False, False, True, 2,
                                      _result_path + '/test_output_app_analysis', 'pdf')
        pydairapp = PyDAIRAPPStats(pydair_args)
        pydairapp.write_freq()
        pydairapp.write_cdr3_len_freq()
        pydairapp.write_diversity_results('vdj')
        pydairapp.stats.get_rarefaction_result()
        pydairapp.plot_figures()
        pydairapp.create_report()
        
    
if __name__ == '__main__':
    unittest.main()

