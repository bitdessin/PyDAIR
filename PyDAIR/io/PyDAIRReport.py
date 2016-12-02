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
    
    def __init__(self, stats, file_path):
        """PyDAIRReport class initialize method.
        
        Args:
            stats (PyDAIRStats): A PyDAIRStats class object.
        """
        
        self.__tmpl_path = str(os.path.dirname(__file__)) + '/../templates'
        self.__env   = Environment(loader = FileSystemLoader(self.__tmpl_path))
        self.__tmpl  = self.__env.get_template('report.html')
        self.__stats = stats
        self.__file_path = file_path
        
    
    
        
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
        freq_stats = [{'name': 'V', 'freq': []},
                       {'name': 'D', 'freq': []},
                       {'name': 'J', 'freq': []}]
        cdr3_stats = []
        rare_stats = []
        vdjidx = {'v': 0, 'd': 1, 'j': 2}
        for i in range(len(self.__stats.samples)):
            bsample = self.__stats.samples.get_record(i)
            sample_stats.append({'name': bsample.name,
                                 'libsize': bsample.len(),
                                 'vfreq': self.__file_path['vfreq'][i],
                                 'dfreq': self.__file_path['dfreq'][i],
                                 'jfreq': self.__file_path['jfreq'][i],
                                 'vdjfreq': self.__file_path['vdjfreq'][i],
                                 'vdjrarefaction': self.__file_path['vdjrarefaction'][i],
                                 'cdr3protlen': self.__file_path['cdr3protlen'][i],
                                 'cdr3nucllen': self.__file_path['cdr3nucllen'][i]})
            
            for g in ['v', 'd', 'j']:
                freq_g = bsample.get_freq(g, prob = True)
                freq_stats[vdjidx[g]]['freq'].append({'sample_name': '"' + bsample.name + '"',
                         'gene': '["' + '","'.join(map(str, freq_g.index.get_values())) + '"]',
                         'freq': '[' + ','.join(map(str, freq_g.tolist())) + ']'})
            
            cdr3_g = bsample.get_freq('cdr3_prot_len', prob = True, remove_cdr3_zero_len = True)
            cdr3_stats.append({'sample_name': '"' + bsample.name + '"',
                                'length': '[' + ','.join(map(str, cdr3_g.index.get_values())) + ']',
                                'freq': '[' + ','.join(map(str, cdr3_g.tolist())) +  ']'})
            
            if bsample.div.rarefaction['vdj'] is not None:
                rare_stats.append({'sample_name': '"' + bsample.name + '"',
                                   'vdj_combn': '[' + ','.join(map(str, bsample.div.rarefaction['vdj'].mean(axis = 1))) + ']',
                                   'sample_size': '[' + ','.join(map(str, bsample.div.rarefaction['vdj'].index.get_values())) + ']'})
            
        report_data['sample_stats'] = sample_stats
        report_data['freq_stats'] = freq_stats
        report_data['cdr3_stats'] = cdr3_stats
        report_data['rare_stats'] = rare_stats
        
        


        
        
        html = self.__tmpl.render(report_data)
        
        report_html = open(file_path, 'w')
        report_html.write(html.encode('utf-8'))
        report_html.close()
        
        
    
    
