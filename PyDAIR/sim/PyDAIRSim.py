import re
import os
import sys
import math
import random
import numpy as np
import pandas as pd
from Bio import SeqIO
from PyDAIR.io.PyDAIRIO import *
from PyDAIR.utils.PyDAIRUtils import *


class PyDAIRSimGeneSet:
    """PyDAIRSimGeneSet class.
    
    The class is for sampling sequences of single of V, D or J gene.
    The sequence name and sequence which are read from FASTA file
    are saved in this class. In addition, the sampling parameters are
    saved in this class.
    
    """
    
    def __init__(self, seq_name, seq_seq, prob, d5_len, d3_len):
        self.name = seq_name
        self.seq = seq_seq
        self.prob = prob
        self.d5_len = d5_len
        self.d3_len = d3_len

    def generate_seqs(self, n, seed = None):
        """Generate V/D/J segments.
        
        Args:
            n (int): The number of sequences should be generated.
            seed (int): The seed used for random sampling.
        """
        if seed is not None:
            np.random.seed(seed)
        
        sampled_id = list(pd.Series(self.name).sample(n, replace = True, weights = self.prob).index)
        del5len = np.random.poisson(self.d5_len, n)
        del3len = np.random.poisson(self.d3_len, n)
        
        sampled_name = []
        sampled_seq = []
        sampled_seq_5del = []
        sampled_seq_3del = []
        for i in range(len(sampled_id)):
            cut_s = del5len[i]
            cut_e = len(self.seq[sampled_id[i]]) - del3len[i]
            sampled_name.append(self.name[sampled_id[i]])
            sampled_seq.append(self.seq[sampled_id[i]][cut_s:cut_e])
            sampled_seq_5del.append(self.seq[sampled_id[i]][:cut_s])
            sampled_seq_3del.append(self.seq[sampled_id[i]][cut_e:])
        return [sampled_name, sampled_seq, sampled_seq_5del, sampled_seq_3del]
        
        

class PyDAIRSimInsSet:
    """PyDAIRSimInsSet class.
    
    The class for saving the parameters to create untemplated nucleotides.
    
    """
    
    def __init__(self, ins_len,
                 sample_set = ['A', 'C', 'G', 'T'],
                 sample_prob = None):
        if sample_prob:
            sample_prob = pd.Series([float(1) / len(sample_set)] * len(sample_set))
        self.ins_len = ins_len
        self.sample_set = pd.Series(sample_set)
        self.prob = sample_prob


    def generate_seqs(self, n, seed = None):
        """Generate untemplated sequences.
        
        Args:
            n (int): The number of sequences should be generated.
            seed (int): The seed used for random sampling.
        
        """
        if seed is not None:
            np.random.seed(seed)
        
        sampled_length = np.random.poisson(self.ins_len, n)
        nontemplate_nucl = []
        
        for _l in sampled_length:
            nontemplate_nucl.append(self.sample_set.sample(_l, replace = True, weights = self.prob).str.cat(sep = ''))
        
        return nontemplate_nucl





class PyDAIRSim:
    """PyDAIRSim class.
    
    """
    
    def __init__(self, v, d, j, vd_ins, dj_ins, p_mutation):
        """Set up paramters for generation of artificial IgH sequences.
        
        """
        self.v = v
        self.d = d
        self.j = j
        self.vd_ins = vd_ins
        self.dj_ins = dj_ins
        self.p_mutation = p_mutation
    
    
    def generate_seqs(self, n = 10000, seed = None, file_name = None):
        """Generate of artificial IgH sequences.
        
        """
        
        if file_name is None:
            raise ValueError('Need \'file_name\'.')
        
        
        # create V/D/J segments
        v_name, v_seq, v_seq_5del, v_seq_3del = self.v.generate_seqs(n, seed)
        d_name, d_seq, d_seq_5del, d_seq_3del = self.d.generate_seqs(n, seed)
        j_name, j_seq, j_seq_5del, j_seq_3del  = self.j.generate_seqs(n, seed)
        
        # create VD, DJ inserted segments
        vd_ins_seqs =  self.vd_ins.generate_seqs(n, seed)
        dj_ins_seqs =  self.dj_ins.generate_seqs(n, seed)
        
        
        with open(file_name, 'w') as fh:
            for i in range(n):
                header = 'SEQ-%09d' % (i + 1) + '|' +  \
                          v_name[i] + '|' + d_name[i] + '|' + j_name[i] + \
                         '|V5DEL:' + v_seq_5del[i] + '|V3DEL:' + v_seq_3del[i] + \
                         '|D5DEL:' + d_seq_5del[i] + '|D3DEL:' + d_seq_3del[i] + \
                         '|J5DEL:' + j_seq_5del[i] + '|J3DEL:' + j_seq_3del[i] + \
                         '|VDINS:' + vd_ins_seqs[i] + '|DJINS:' + dj_ins_seqs[i]
                seq = v_seq[i] + vd_ins_seqs[i] + d_seq[i] + dj_ins_seqs[i] + j_seq[i]
                # add mutation into each sequence
                seq = self.__mutate_seq(seq, self.p_mutation)
                fh.write('>' + header + '\n' + seq + '\n')
    
    
    def __mutate_seq(self, seq, prob, seed = None):
        if seed is not None:
            np.random.seed(seed)
        
        b = ['A', 'C', 'G', 'T', ' ']
        probs = np.random.rand(len(seq))
        
        for i in range(len(probs)):
            if probs[i] < prob:
                seq = seq[:i] + np.random.choice(b) + seq[(i+1):]
        seq = seq.replace(' ', '')
        return seq



class PyDAIRSimEval:
    """PyDAIRSimEval class.
    
    """
    
    def __init__(self, sim_file, pydair_file):
        """Set up parameters for evaluating performances.
        
        Args:
            sim_file (str): A file path to FASTA file generated by PyDAIR sim mode.
            pydair_file (str): A file path to PYDAIR file generated by PyDAIR parse mode.
        
        """
        
        self.sim_file = sim_file
        self.pydair_file = pydair_file
    
    
    def eval(self, file_name):
        """Calculate performances.
        
        Args:
            file_name (str): A file path to save the evaluated results.
        
        """
        
        utilsobj = PyDAIRUtils()
        
        # loading the simulation conditions
        seq_id = []
        trueobj = {}
        with open(self.sim_file, 'r') as fa_fh:
            for record in SeqIO.parse(fa_fh, 'fasta'):
                _seq_id, _v, _d, _j, \
                _v5del, _v3del, _d5del, _d3del, _j5del, _j3del, \
                _vdins, _djins = record.id.split('|')
                
                seq_id.append(record.id)
                
                trueobj[record.id] = {
                    'v': _v, 'd': _d, 'j': _j,
                    'v5del': _v5del.replace('V5DEL:', ''),
                    'v3del': _v3del.replace('V3DEL:', ''),
                    'd5del': _d5del.replace('D5DEL:', ''),
                    'd3del': _d3del.replace('D3DEL:', ''),
                    'j5del': _j5del.replace('J5DEL:', ''),
                    'j3del': _j3del.replace('J3DEL:', ''),
                    'vdins': _vdins.replace('VDINS:', ''),
                    'djins': _djins.replace('DJINS:', '')
                }
        
        # loading the analysis results.
        estobj = {}
        pydairfh = PyDAIRIO(self.pydair_file, 'r')
        for igseq in pydairfh.parse():
            estobj[igseq.query.name] = {
                'v': utilsobj.none_to_dot(igseq.v.sbjct.name),
                'd': utilsobj.none_to_dot(igseq.d.sbjct.name),
                'j': utilsobj.none_to_dot(igseq.j.sbjct.name),
                'v3del': utilsobj.none_to_dot(igseq.indels.v_deletion),
                'j5del': utilsobj.none_to_dot(igseq.indels.j_deletion),
                'vjins': utilsobj.none_to_dot(igseq.indels.vj_insertion)
            }
        
        # calculate performances
        v_results = []
        d_results = []
        j_results = []
        v3del_results = []
        j5del_results = []
        
        outfh = open(file_name, 'w')
        for _seq_id in seq_id:
            if _seq_id not in estobj:
                estobj[_seq_id] = {
                    'v': '.', 'd': '.', 'j': '.',
                    'v3del': '.', 'j5del': '.', 'vjins': '.'
                }
            _v, _vt = self.__eval_est_result(trueobj, estobj, _seq_id, 'v')
            _d, _dt = self.__eval_est_result(trueobj, estobj, _seq_id, 'd')
            _j, _jt = self.__eval_est_result(trueobj, estobj, _seq_id, 'j')
            _v3del, _v3delt = self.__eval_est_result(trueobj, estobj, _seq_id, 'v3del')
            _j5del, _j5delt = self.__eval_est_result(trueobj, estobj, _seq_id, 'j5del')
            
            outfh.write('\t'.join([_seq_id.split('|')[0],
                                   _v, _d, _j, _v3del, _j5del]) + '\n')
            
            v_results.append(_vt)
            d_results.append(_dt)
            j_results.append(_jt)
            v3del_results.append(_v3delt)
            j5del_results.append(_j5delt)
        
        sumtxt = ['Feature\tCorrect\tIncorrect\tUnidentifiable\tTotal']
        sumtxt.append(self.__eval_calc_sum('V', v_results))
        sumtxt.append(self.__eval_calc_sum('D', d_results))
        sumtxt.append(self.__eval_calc_sum('J', j_results))
        sumtxt.append(self.__eval_calc_sum('V3DEL', v3del_results))
        sumtxt.append(self.__eval_calc_sum('J5DEL', j5del_results))
        
        outfh.write('\n\n'+ '\n'.join(sumtxt) + '\n')
        
        outfh.close()
            
            
    def __eval_calc_sum(self, g, x):
        y = pd.Series(x)
        s = [g, str(sum(y == 'T')), str(sum(y == 'F')),
                str(sum(y == 'U')), str(y.size)]
        return('\t'.join(s))
    
        
    def __eval_est_result(self, trueobj, estobj, seq_id, k):
        tag = self.__get_eval_tag(trueobj[seq_id][k], estobj[seq_id][k])
        evaltxt = trueobj[seq_id][k] + '\t' + estobj[seq_id][k] + '\t' +  tag
        return [evaltxt, tag]
        
       
     
    def __get_eval_tag(self, true, est):
        if str(est) == '.':
            return 'U'
        if str(true) == str(est):
            return 'T'
        else:
            return 'F'
        






