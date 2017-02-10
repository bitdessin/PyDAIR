import os
import unittest
from PyDAIR.sim.PyDAIRSim import *
from PyDAIR.app.PyDAIRAPP import *
from PyDAIR.utils.PyDAIRUtils import *
from Bio import SeqIO

_data_path = os.path.join(os.path.dirname(__file__), 'data/samples')
_db_path = os.path.join(os.path.dirname(__file__), 'data/db')
_result_path = os.path.join(os.path.dirname(__file__), 'data/results')



class Test_pydair_simulation(unittest.TestCase):
    def setUp(self):
        self.n = 100
        self.v_fa = _db_path + '/v.fa'
        self.d_fa = _db_path + '/d.fa'
        self.j_fa = _db_path + '/j.fa'
        self.v_5del = 10
        self.v_3del = 3
        self.d_5del = 3
        self.d_3del = 3
        self.j_5del = 5
        self.j_3del = 10
        self.mutate = 0.05
        self.vd_ins = 5
        self.dj_ins = 5
        self.seed = 1010
        self.o_file_name = _result_path + '/test_generatesimseq.fa'
        self.i_file_name = _data_path + '/sample.simseq.fa'
    
    
    def test_pydair_sim(self):
        # read V gene fasta
        v_name = []
        v_seq = []
        with open(self.v_fa, 'r') as fa_fh:
            for record in SeqIO.parse(fa_fh, 'fasta'):
                v_name.append(record.id)
                v_seq.append(str(record.seq).upper())
        # read J gene fasta
        j_name = []
        j_seq = []
        with open(self.j_fa, 'r') as fa_fh:
            for record in SeqIO.parse(fa_fh, 'fasta'):
                j_name.append(record.id)
                j_seq.append(str(record.seq).upper())
        # read D gene fasta
        d_name = []
        d_seq = []
        with open(self.d_fa, 'r') as fa_fh:
            for record in SeqIO.parse(fa_fh, 'fasta'):
                d_name.append(record.id)
                d_seq.append(str(record.seq).upper())
        
        # create object
        v_obj = PyDAIRSimGeneSet(v_name, v_seq, None,
                                 self.v_5del, self.v_3del)
        d_obj = PyDAIRSimGeneSet(d_name, d_seq, None,
                                 self.d_5del, self.d_3del)
        j_obj = PyDAIRSimGeneSet(j_name, j_seq, None,
                                 self.j_5del, self.j_3del)
        vd_obj = PyDAIRSimInsSet(self.vd_ins)
        dj_obj = PyDAIRSimInsSet(self.dj_ins)
        
        simobj = PyDAIRSim(v = v_obj, d = d_obj, j = j_obj,
                           vd_ins = vd_obj, dj_ins = dj_obj,
                           p_mutation = self.mutate)
        
        # create artificial IgH sequences
        simobj.generate_seqs(self.n, self.seed, self.o_file_name)
        
    
    def test_pydair_simeval(self):
        simobj = PyDAIRSimEval(_data_path + '/sample.simseq.fa',
                               _data_path + '/sample.simseq.parsed.pydair')
        simobj.eval(file_name = _result_path + '/test_sim.eval.txt')







if __name__ == '__main__':
    unittest.main()


