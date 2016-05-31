import os
import math
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from PyDAIR.seq.IgSeq import IgSeq
from PyDAIR.io.PyDAIRIO import *

_data_path = os.path.join(os.path.dirname(__file__), 'data')


class Test_pydair_io(unittest.TestCase):
    
    
    def setUp(self):
        # Inputs
        self.pydair_input_path       = _data_path + '/sample.pydair'
        # Outputs
        self.pydair_output_path       = _data_path + '/test_output_io.pydair'
        self.pydairsimple_output_path = _data_path + '/test_output_io.pydair.simple'
    
    
    def test_pydair_io(self):
        pydair_i_fh = PyDAIRIO(self.pydair_input_path, 'r', 'pydair')
        pydair_o_fh = PyDAIRIO(self.pydair_output_path, 'w', 'pydair')
        # read record in PyDAIR flat file
        for igseq in pydair_i_fh.parse():
            # print contents
            print('-----------------------------------------------')
            if igseq.query.orf is not None:
                seq_nucl  = igseq.query.seq[igseq.query.orf:int(math.floor((len(igseq.query.seq) - igseq.query.orf) / 3) * 3 + igseq.query.orf)]
                cdr3_nucl = igseq.query.seq[igseq.variable_region.cdr3[0]:igseq.variable_region.cdr3[1]]
                print(len(seq_nucl))
                seq_prot  = str(Seq(seq_nucl, generic_dna).translate())
                cdr3_prot = str(Seq(cdr3_nucl, generic_dna).translate())
                print('nucl:  ' + seq_prot)
                print('prot:  ' + cdr3_prot)
            else:
                print('nucl:  .')
                print('prot:  .')
            # write object into new file
            pydair_o_fh.write(igseq)
        print('-----------------------------------------------')
        pydair_i_fh.close()
        pydair_o_fh.close()
    
    
    def test_pydairsimple_io(self):
        pydair_i_fh = PyDAIRIO(self.pydair_input_path, 'r', 'pydair')
        pydair_o_fh = PyDAIRIO(self.pydairsimple_output_path, 'w', 'simple')
        for igseq in pydair_i_fh.parse():
            pydair_o_fh.write(igseq)
        pydair_i_fh.close()
        pydair_o_fh.close()
    
    
if __name__ == '__main__':
    unittest.main()

