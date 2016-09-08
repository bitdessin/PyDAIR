from __future__ import division
import re
import os
import sys
import math
import random
import numpy as np
import pandas as pd
from PyDAIR.seq.IgSeq import *
from PyDAIR.io.PyDAIRIO import *


logging.basicConfig(level = logging.INFO, format = '%(levelname)-8s %(message)s')


class PyDAIRDiversity:
    def __init__(self):
        self.rarefaction        = {'vdj': None, 'cdr3': None}
        self.samplingresampling = {'vdj': None, 'cdr3': None}




class PyDAIRStatsRecord:
    '''
    The class for storing the summarised PyDAIR data.
    
    One PyDAIR file should has one BLGIHStatsRecord class object. If there
    are more than one PyDAIR files, use a number of PyDAIRStatsRecord class
    objects to save these data, and save all PyDAIRStatsRecord class objects
    into PyDAIRStatsRecords class.
    This class requires the list of V gene names, D gene names, J gene names,
    CDR3 nucleotide and protein seuqence, and the stop codon tags.
    '''
    def __init__(self, name = None, v = None, d = None, j = None,
                 cdr3_nucl_seq = None, cdr3_prot_seq = None, stop_codon_tag = None,
                 contain_ambiguous_D = False, contain_stopcodon = False):
        # calculate the number of entries
        
        # set default data
        self.name = name
        self.vdj  = pd.DataFrame({'v': v, 'd': d, 'j': j}, columns = ['v', 'd', 'j']).fillna(value = np.nan)
        self.cdr3 = pd.DataFrame({'nucl_seq': cdr3_nucl_seq, 'prot_seq': cdr3_prot_seq,
                                  'nucl_len': pd.Series(cdr3_nucl_seq).str.len(),
                                  'prot_len': pd.Series(cdr3_prot_seq).str.len()},
                                  columns = ['nucl_seq', 'prot_seq', 'nucl_len', 'prot_len']).fillna(value = np.nan)
        self.stop_codon_tag = pd.Series(stop_codon_tag).fillna(value = np.nan)
        
        # set filters
        filter_ambigoD = self.vdj.d.notnull()
        filter_stcodon = pd.Series(self.stop_codon_tag == 'N')
        
        if contain_ambiguous_D:
            filter_ambigoD = pd.Series([True] * self.vdj.shape[0])
        if contain_stopcodon:
            filter_stcodon = pd.Series([True] * self.cdr3.shape[0])

        filters = pd.Series(filter_ambigoD & filter_stcodon)
        
        # filter data
        self.vdj  = self.vdj[filters]
        self.cdr3 = self.cdr3[filters]
        self.stop_codon_tag = self.stop_codon_tag[filters]
        
        # diversity study
        self.div = PyDAIRDiversity()
    
    
    
    def __len__(self):
        return self.vdj.shape[0]
    
    
    def len(self):
        return self.__len__
    
    
    def get_freq(self, gene):
        freq = None
        if gene.lower() == 'cdr3_prot_len':
            freq = self.cdr3.prot_len.value_counts(dropna = False)
            freq = freq.sort_index(ascending = True)
            freq.name = 'frequency'
        elif gene.lower() == 'cdr3_nucl_len':
            freq = self.cdr3.nucl_len.value_counts(dropna = False)
            freq = freq.sort_index(ascending = True)
            freq.name = 'frequency'
        elif gene.lower() == 'v':
            freq = self.vdj.v.value_counts(dropna = False)
            freq.name = 'frequency'
        elif gene.lower() == 'd':
            freq = self.vdj.d.value_counts(dropna = False)
            freq.name = 'frequency'
            freq.index = ['unidentifiable' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
        elif gene.lower() == 'j':
            freq = self.vdj.j.value_counts(dropna = False)
            freq.name = 'frequency'
        elif gene.lower() == 'vdj':
            __sep = '__________'
            vdj_combinations = self.vdj.v.replace(np.nan, 'NaN') + __sep + \
                               self.vdj.d.replace(np.nan, 'NaN') + __sep + \
                               self.vdj.j.replace(np.nan, 'NaN')
            freq = vdj_combinations.value_counts(dropna = False)
            vdj_v = []
            vdj_d = []
            vdj_j = []
            for vdj_combination in freq.index:
                vv, dd, jj = vdj_combination.split(__sep)
                vdj_v.append(vv)
                vdj_d.append(dd)
                vdj_j.append(jj)
            freq = pd.DataFrame({'v': vdj_v, 'd': vdj_d, 'j': vdj_j, 'frequency': freq.values},
                                columns = ['v', 'd', 'j', 'frequency']).replace('NaN', np.nan)
        else:
            raise ValueError('The \'gene\' argument of \'get_freq\' gene should be one of \'v\', \'d\', \'j\', and \'vdj\'.')
        
        return freq
    
    
    def samplingresampling_study(self, data = None, n = 1000):
        if data is None or data == 'all':
            data = ['vdj', 'cdr3']
        if not isinstance(data, list):
            data = [data]

        for data_i in data:
            if data_i == 'vdj':
                dat = self.vdj
                __sep = '__________'
                dat_for_study = dat.v.replace(np.nan, 'NaN') + __sep + \
                                dat.d.replace(np.nan, 'NaN') + __sep + \
                                dat.j.replace(np.nan, 'NaN')
                dat_for_study = pd.Series(dat_for_study)
                
            if data_i == 'cdr3':
                dat_for_study = pd.Series(self.cdr3.prot_seq)
            
            self.div.samplingresampling[data_i] = self.__samplingresampling_study(dat_for_study, n)
        
        
    def __samplingresampling_study(self, dat, n):
        population       = dat
        population_size  = len(dat)
        resampling_size_ratio = range(101)
        
        # get resampling sizes
        resampling_sizes = []
        for r in resampling_size_ratio:
            resampling_sizes.append(int(round(0.01 * r * population_size)))
        x = [0]   # total sampling sizes
        y = []    # the number of new sampled CDR3 sequence
        for i in range(n):
            y_try_n      = []
            population_n = list(population)
            sampled_n    = set([])
            
            for s in range(len(resampling_sizes) - 1):
                samplingsize = resampling_sizes[s + 1] - resampling_sizes[s]
                if len(population_n) < samplingsize * 1.5:
                    samplingsize = len(population_n)
                if len(population_n) <= samplingsize:
                    sampled_idx = range(0, len(population_n))
                else:
                    sampled_idx = random.sample(range(0, len(population_n)), samplingsize)
                s_left = []
                s_smpl = []
                for smplidx_in_population_n in range(len(population_n)):
                    if smplidx_in_population_n in sampled_idx:
                        s_smpl.append(population_n[smplidx_in_population_n])
                    else:
                        s_left.append(population_n[smplidx_in_population_n])
                sampled_items = s_smpl
                population_n  = s_left
                sampled_items = set(sampled_items)
                y_try_n.append(len(sampled_items.difference(sampled_n)))
                sampled_n = sampled_n.union(sampled_items)
                if i == 0:
                    x.append(x[s] + samplingsize)
            y.append(y_try_n)
        x.pop(0)
        y = pd.DataFrame(y, columns = resampling_size_ratio[1:]).T
        y.columns = ['try_' + str(n_try + 1) for n_try in range(n)]
        return y
    
    
    
    def rarefaction_study(self, data = None, n = 1000):
        if data is None or data == 'all':
            data = ['vdj', 'cdr3']
        if not isinstance(data, list):
            data = [data]
        
        for data_i in data:
            if data_i == 'vdj':
                dat = self.vdj
                __sep = '__________'
                dat_for_study = dat.v.replace(np.nan, 'NaN') + __sep + \
                                dat.d.replace(np.nan, 'NaN') + __sep + \
                                dat.j.replace(np.nan, 'NaN')
                dat_for_study = pd.Series(dat_for_study)
                
            if data_i == 'cdr3':
                dat_for_study = pd.Series(self.cdr3.prot_seq)
                
            self.div.rarefaction[data_i] = self.__rarefaction_study(dat_for_study, n)
        
        
    

    def __rarefaction_study(self, dat, n):
        population_size      = len(dat)
        sampling_sizes       = self.__get_sampling_sizes(population_size)
        sampled_unique_items = [None] * n
        
        # rarefaction study
        for i in range(n):
            sampled_unique_items_i = []
            for s in sampling_sizes:
                sampled_id = random.sample(range(0, population_size), s)
                sampled_unique_items_i.append(dat[sampled_id].unique().shape[0])
            sampled_unique_items[i] = sampled_unique_items_i
        
        # set results into DataFrame class object
        sampled_unique_items = pd.DataFrame(sampled_unique_items, columns = sampling_sizes).T
        sampled_unique_items.columns = ['try_' + str(n_try + 1) for n_try in range(n)]
        
        return sampled_unique_items
    
    

    def __get_sampling_sizes(self, sample_size):
        sampling_sizes = []
        n_digits = int(math.log10(sample_size) + 1)
        # 1st sampling sizes
        for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            sampling_sizes.append(i * 10 ** (n_digits - 1 - 2))
        # 2nd sampling sizes
        for j in [2, 4, 6, 8, 10]:
            sampling_sizes.append(j * 10 ** (n_digits - 1 - 1))
        # 3rd sampling sizes
        k = 2
        while max(sampling_sizes) < sample_size:
            if k * 10 ** (n_digits - 1) < sample_size:
                sampling_sizes.append(k * 10 ** (n_digits - 1))
            else:
                sampling_sizes.append(sample_size)
            k += 1
        return sampling_sizes


       
        
   
class PyDAIRStatsRecords:
    '''PyDAIRStatsRecords for storing the PyDAIR data.
    
    This class provides the list of PyDAIRStatsRecord and some function of list.
    '''
    def __init__(self):
        self.__records = []
        self.__index = 0
    
    def __len__(self):
        return len(self.__records)
    
    def len(self):
        return self.__len__()

    def __iter__(self):
        self.__index = 0
        return self
    
    def next(self):
        try:
            stats_record = self.__records[self.__index]
        except IndexError:
            raise StopIteration
        self.__index += 1
        return stats_record
    
    def append(self, stats_record):
        self.__records.append(stats_record)
    
    def get_record(self, i):
        try:
           stats_record = self.__records[i]
        except IndexError:
            raise IndexError
        return stats_record
    
    def set_record(self, i, stats_record):
        self.__records[i] = stats_record


    
class PyDAIRStats:
    '''
    Statistics analysis of repertoire sequences.
    
    The class for storing the PyDAIR data file. Before use initialize this class,
    one should run create PyDAIR format files by other PyDAIR functions. This class
    stored many data of analyzed data.
    '''
    def __init__(self, pydair_file, pydair_format = None, pydair_id = None,
                 contain_ambiguous_D = False, contain_stopcodon = False):
        if pydair_format is None:
            pydair_format = 'PyDAIR'
        if pydair_id is None:
            pydair_id = []
            for i in range(len(pydair_file)):
                pydair_id.append('individual ' + str(i + 1))
        
        self.__pydair_file   = pydair_file
        self.__pydair_format = pydair_format
        self.__pydair_id     = pydair_id
        self.__contain_ambiguous_D = contain_ambiguous_D
        self.__contain_stopcodon = contain_stopcodon
        self.samples   = None
        
        # parse PyDAIR files
        self.__parse_pydair_files()
        
    
    def __parse_pydair_files(self):
        self.samples  = PyDAIRStatsRecords()
        for i in range(len(self.__pydair_file)):
            v = []
            d = []
            j = []
            cdr3_prot_seq = []
            cdr3_nucl_seq = []
            stop_codon_tag = []
            
            pydair_fh = PyDAIRIO(self.__pydair_file[i], 'r', self.__pydair_format)
            for igseq in pydair_fh.parse():
                v.append(igseq.v.sbjct.name)
                d.append(igseq.d.sbjct.name)
                j.append(igseq.j.sbjct.name)
                
                cdr3_data = igseq.get_cdr3_data()
                cdr3_prot_seq.append(cdr3_data.prot_seq)
                cdr3_nucl_seq.append(cdr3_data.nucl_seq)
                stop_codon_tag.append(cdr3_data.stop_codon_tag)
            
            sample_record = PyDAIRStatsRecord(self.__pydair_id[i], v, d, j,
                                              cdr3_nucl_seq, cdr3_prot_seq, stop_codon_tag,
                                              self.__contain_ambiguous_D, self.__contain_stopcodon)
            self.samples.append(sample_record)
    
    
    
    
    
