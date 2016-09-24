import math
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from PyDAIR.seq.IgSeq import IgSeq
from PyDAIR.io.PyDAIRIO import *
from PyDAIR.stats.PyDAIRStats import *
from PyDAIR.plot.PyDAIRPlot import *

_data_path = os.path.join(os.path.dirname(__file__), 'data')
_result_path = os.path.join(os.path.dirname(__file__), 'data/results')


class Test_pydair_plot(unittest.TestCase):
    
    
    def setUp(self):
        # input data
        pydair_files = [_data_path + '/sample.1.pydair',
                        _data_path + '/sample.2.pydair',
                        _data_path + '/sample.3.pydair']
        pydair_id    = ['Sample 1', 'Sample 2', 'Sample 3']
        
        # output settings
        self.freq_v_png = _result_path + '/test_output_plot.freq.v.png'
        self.freq_d_png = _result_path + '/test_output_plot.freq.d.png'
        self.freq_j_png = _result_path + '/test_output_plot.freq.j.png'
        self.prob_v_pdf = _result_path + '/test_output_plot.prob.v.pdf'
        self.prob_d_pdf = _result_path + '/test_output_plot.prob.d.pdf'
        self.prob_j_pdf = _result_path + '/test_output_plot.prob.j.pdf'
        
        self.dist_cdr3_png_1 = _result_path + '/test_output_plot.cdr3len.1.png'
        self.dist_cdr3_png_2 = _result_path + '/test_output_plot.cdr3len.2.png'
        self.dist_cdr3_png_3 = _result_path + '/test_output_plot.cdr3len.3.png'
        
        # create stat object
        self.stats   = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        self.stats_D = PyDAIRStats(pydair_files, 'pydair', pydair_id, contain_ambiguous_D = True)
        
    
    def test_plot_gene_freq(self):
        bplots   = PyDAIRPlot(self.stats, 'ggplot')
        bplots.barplot_freq(gene = 'v', fig_name = self.freq_v_png, fig_format = 'png')
        bplots.barplot_freq(gene = 'd', fig_name = self.freq_d_png, fig_format = 'png',
                            gene_names = ['Dm1', 'Dm2', 'Dm3', 'Dm4', 'Dm5', 'Dm6', 'Dm7'])
        bplots.barplot_freq(gene = 'j', fig_name = self.freq_j_png, fig_format = 'png',
                            main = 'J gene usage frequency', xlab = 'Gene', ylab = 'Counts')
        
        bplots_D = PyDAIRPlot(self.stats_D, 'fivethirtyeight')
        bplots_D.barplot_freq(gene = 'v', fig_name = self.prob_v_pdf, fig_format = 'pdf', prob = True)
        bplots_D.barplot_freq(gene = 'd', fig_name = self.prob_d_pdf, fig_format = 'pdf', prob = True)
        bplots_D.barplot_freq(gene = 'j', fig_name = self.prob_j_pdf, fig_format = 'pdf', prob = True)
    
    
    
    def test_plot_cdr3len(self):
        bplots   = PyDAIRPlot(self.stats, 'classic')
        bplots.hist_cdr3_len(fig_name = self.dist_cdr3_png_2,
                             xlim = [10, 30], main = 'CDR3 AA Length')
        bplots.hist_cdr3_len(fig_name = self.dist_cdr3_png_3,
                             main = 'CDR3 Length Distributions', prob = True)
        
        bplots_D = PyDAIRPlot(self.stats_D, 'bmh')
        bplots_D.hist_cdr3_len(fig_name = self.dist_cdr3_png_1, xlim = [5, 25])



if __name__ == '__main__':
    unittest.main()




