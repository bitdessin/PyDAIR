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
            file_path (dictionary): A dictionary to the summary files.
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
        for i in range(len(self.__stats.samples)):
            bsample = self.__stats.samples.get_record(i)
            _s = {'name': bsample.name, 'libsize': bsample.len()}
            _s.update(self.__file_path[bsample.name])
            sample_stats.append(_s)
        _s = {'name': 'all', 'libsize': '-'}
        _s.update(self.__file_path['all'])
        sample_stats.append(_s)
        report_data['sample_stats'] = sample_stats
        
        
        freq_stats = []
        for g in ['V', 'D', 'J']:
            _s = {'freq_name': g, 'freq_title': g, 'freq': []}
            for i in range(len(self.__stats.samples)):
                _b = self.__stats.samples.get_record(i)
                _f = _b.get_summary(g.lower(), prob = True)
                _s['freq'].append({
                    'sample_name': '"' + _b.name + '"',
                    'gene': '["' + '","'.join(map(str, _f.index.get_values())) + '"]',
                    'freq': '[' + ','.join(map(str, _f.tolist())) + ']'
                })
            freq_stats.append(_s)
        report_data['freq_stats'] = freq_stats
        
        
        dist_stats = []
        dist_stats_names = ['cdr3_prot_len', 'v_del_len', 'j_del_len', 'vj_ins_len']
        dist_stats_title = {'cdr3_prot_len': 'CDR3 amino acid length',
                            'v_del_len': '5\'-end V deletion length',
                            'j_del_len': '3\'-end J deletion length',
                            'vj_ins_len': 'V-J junction insertion length'}
        for g in dist_stats_names:
            _s = {'dist_name': g, 'dist_title': dist_stats_title[g], 'freq': []}
            for i in range(len(self.__stats.samples)):
                _b = self.__stats.samples.get_record(i)
                _f = _b.get_summary(g, prob = True)
                _s['freq'].append({
                    'sample_name': '"' + _b.name + '"',
                    'len': '["' + '","'.join(map(str, _f.index.get_values())) + '"]',
                    'freq': '[' + ','.join(map(str, _f.tolist())) + ']'
                })
            dist_stats.append(_s)
        report_data['dist_stats'] = dist_stats
        
        
        
        
        rare_stats = []
        vdjidx = {'v': 0, 'd': 1, 'j': 2}
        for i in range(len(self.__stats.samples)):
            bsample = self.__stats.samples.get_record(i)
            if bsample.div.rarefaction['vdj'] is not None:
                rare_stats.append({'sample_name': '"' + bsample.name + '"',
                                   'vdj_combn': '[' + ','.join(map(str, bsample.div.rarefaction['vdj'].mean(axis = 1))) + ']',
                                   'sample_size': '[' + ','.join(map(str, bsample.div.rarefaction['vdj'].index.get_values())) + ']'})
            
        report_data['rare_stats'] = rare_stats
        
        
        
        
        html = self.__tmpl.render(report_data)
        
        report_html = open(file_path, 'w')
        report_html.write(html.encode('utf-8'))
        report_html.close()
        
        
    
    
