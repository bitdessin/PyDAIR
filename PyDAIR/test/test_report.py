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

_data_path = os.path.join(os.path.dirname(__file__), 'data')


class Test_bin(unittest.TestCase):
    
    def setUp(self):
        
        # set up statistic data
        pydair_files = [_data_path + '/sample.1.pydair',
                        _data_path + '/sample.2.pydair',
                        _data_path + '/sample.3.pydair']
        pydair_id    = ['Sample 1', 'Sample 2', 'Sample 3']
        self.stats   = PyDAIRStats(pydair_files, 'pydair', pydair_id)
        
        # set up template
        self.__report_1 = _result_path + '/test_report_1.html'
        self.__report_2 = _result_path + '/test_report_2.html'
        
        
    def test_report_1(self):
        stats = self.stats
        stats.rarefaction_study('vdj', n = 10)
        report = PyDAIRReport(stats)
        report.render(self.__report_1)
    
    def test_report_2(self):
        stats = self.stats
        report = PyDAIRReport(stats)
        report.render(self.__report_2)
    
        

if __name__ == '__main__':
    unittest.main()

