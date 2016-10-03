import re
import os
import sys
import math
from PyDAIR.utils.PyDAIRUtils import *
from PyDAIR.seq.IgSeq import *
from jinja2 import Environment, FileSystemLoader



class PyDAIRReport: 
    """PyDAIR report render class.
    
    
    This class renders the anaysis results that obtained from PyDAIRStats into HTML template.
    It is expected to use PyDAIRReport class methods after analyzing data in PYDAIR format file.
    """
    
    def __init__(self, stats):
        """PyDAIRReport class initialize method.
        
        Args:
            stats (PyDAIRStats): A PyDAIRStats class object.
        """
        
        self.__tmpl_path = str(os.path.dirname(__file__)) + '/../templates'
        self.__env   = Environment(loader = FileSystemLoader(self.__tmpl_path))
        self.__tmpl  = self.__env.get_template('report.html')
        self.__stats = stats
        self.__data_path = None
        
    
    def set_data_path(self, sample, v_path = None, d_path = None, j_path = None,
                            vdj_path = None, cdr3len = None, cdr3fa = None):
        self.__data_path_v = v_path
        self.__data_path_d_freq = None
        self.__data_path_j_freq = None
        self.__data_path_vdj_freq = None
        self.__data_path_vdj_freq = None
        self.__data_path_cdr3len  = None
        self.__data_path_cdr3 = None
        
    
        
    def render(self, file_path):
        """Render data into template.
        
        Args:
            file_path (str): A file name as string to write to.
        
        Render all analysis results into HTML template using Jinja2.
        The template (report.html) is saved in the templates directory of this package.
        """
        
        report_data = {}
        
        # set up report title
        report_data['title'] = 'PyDAIR Analysis Report'
        
        
        # set up basic statistics
        sample_stats = []
        for sample in self.__stats.samples:
            sample_stats.append({'name': sample.name,
                                 'libsize': sample.len()})
        report_data['sample_stats'] = sample_stats
        
        
        # set up V, D, and J statistics
        freq_stats = []
        for gene in ['v', 'd', 'j']:
            freq_g = self.__stats.get_freq(gene)
            freq_stats.append({'gene': gene.upper(),
                               'count': min([freq_g.shape[0], 10]) + 2,
                               'freq_csv' :  freq_g.to_csv(index_label = 'Gene', sep = '\t'),
                               'freq_json':  freq_g.to_json()})
        report_data['freq_stats'] = freq_stats
        
        
        # set up CDR3
        cdr3_len_freq = self.__stats.get_cdr3len_freq()
        cdr3_stats = {'freq_csv': cdr3_len_freq.to_csv(index_label = "Length", sep = '\t'),
                      'freq_json': cdr3_len_freq.to_json()}
        report_data['cdr3_stats'] = cdr3_stats
        
        
        # set up VDJ rarefaction data
        rfdata = self.__stats.get_rarefaction_result(fun = 'mean')
        if rfdata is not None:
            report_data['rarefaction_stats'] = {'tries_csv': rfdata.to_csv(index_label = "CapturedSeq", sep = '\t'),
                                                'tries_json': rfdata.to_json()}
        else:
            report_data['rarefaction_stats'] = None
        
        
        
        html = self.__tmpl.render(report_data)
        
        report_html = open(file_path, 'w')
        report_html.write(html.encode('utf-8'))
        report_html.close()
        
        
    
    
