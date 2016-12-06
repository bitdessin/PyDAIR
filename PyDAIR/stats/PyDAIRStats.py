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




class PyDAIRDiversity:
    """Class for saving diversity study results.
    
    """
    
    def __init__(self):
        """PyDAIRDiversity class initialize method.
        
        Set up ``None`` object in default.
        """
        
        self.rarefaction        = {'vdj': None, 'cdr3': None}
        self.samplingresampling = {'vdj': None, 'cdr3': None}




class PyDAIRStatsRecord:
    """The class for storing the summarised PyDAIR data.
    
    One PyDAIR file should has one BLGIHStatsRecord class object. If there
    are more than one PyDAIR files, use a number of PyDAIRStatsRecord class
    objects to save these data, and save all PyDAIRStatsRecord class objects
    into PyDAIRStatsRecords class.
    This class requires the list of V gene names, D gene names, J gene names,
    CDR3 nucleotide and protein seuqence, and the stop codon tags.
    """
    
    def __init__(self, name = None, v = None, d = None, j = None, orf = None,
                 cdr3_nucl_seq = None, cdr3_prot_seq = None,
                 v_del = None, j_del = None, vj_ins = None,
                 discard_ambiguous_D = False, productive_only = False):
        """PyDAIRStatsRecord class initialize method.
        
        Args:
            name (str): Sample name.
            v (list): Assigned V gene names.
            d (list): Assigned D gene names.
            j (list): Assigned J gene names.
            orf (list): A list of ORF.
            cdr3_nucl_seq (list): A list of CDR3 nucleotide sequences.
            cdr3_prot_seq (list): A list of CDR3 amino acid sequences.
            v_del (list): Deleted nucleotides of 3'-end V gene.
            j_del (list): Deleted nucleotides of 5'-end J gene.
            vj_ins (list): Inserted nucleotides.
            discard_ambiguous_D (bool): If ``True``, ambiguous D gene will be discarded before analysis.
            productive_only (bool): If ``True``, analysis sequence with stop codons.
        """
        # calculate the number of entries
        
        # set default data
        self.name = name
        self.vdj  = pd.DataFrame({'v': v, 'd': d, 'j': j}, columns = ['v', 'd', 'j']).fillna(value = np.nan)
        self.cdr3 = pd.DataFrame({'nucl_seq': cdr3_nucl_seq, 'prot_seq': cdr3_prot_seq,
                                  'nucl_len': pd.Series(cdr3_nucl_seq).str.len(),
                                  'prot_len': pd.Series(cdr3_prot_seq).str.len()},
                                  columns = ['nucl_seq', 'prot_seq', 'nucl_len', 'prot_len']).fillna(value = np.nan)
        self.indels = pd.DataFrame({'v_del': v_del, 'v_del_len': pd.Series(v_del).str.len(),
                                    'j_del': j_del, 'j_del_len': pd.Series(j_del).str.len(),
                                    'vj_ins': vj_ins, 'vj_ins_len': pd.Series(vj_ins).str.len()},
                                    columns = ['v_del', 'v_del_len', 'j_del', 'j_del_len',
                                               'vj_ins', 'vj_ins_len']).fillna(value = np.nan)
        
        
        # set filters
        filter_ambigoD = None
        if discard_ambiguous_D:
            filter_ambigoD = self.vdj.d.notnull()
        else:
            filter_ambigoD = pd.Series([True] * self.vdj.shape[0])
        
        filter_stopcodon = None
        if productive_only:
            filter_stopcodon = pd.Series(orf).notnull()
        else:
            filter_stopcodon = pd.Series([True] * self.vdj.shape[0])
            
        filters = pd.Series(filter_ambigoD & filter_stopcodon)
        
        # filter data
        self.vdj  = self.vdj[filters]
        self.cdr3 = self.cdr3[filters]
        self.indels = self.indels[filters]
        
        # diversity study
        self.div = PyDAIRDiversity()
    
    
    
    def __len__(self):
        return self.vdj.shape[0]
    
    
    def len(self):
        """Retrive the number of sequences analyzed.
        
        """
        
        return self.__len__()
    
    
    def get_summary(self, data_type, prob = False, func = 'mean'):
        """Returns summary statistics.
        
        Args:
            data_type (str): A string to specify data type. ``v``, 
                             ``d``, ``j``, ``vdj``, ``cdr3_nucl_len``,
                             ``cdr3_prot_len``, ``v_del_len``, ``j_del_len``,
                             and ``vj_ins_len`` are supported.
            prob (bool): If ``True``, calculate the probability.
            func (str): A string to specify data calculations when data_type is vdj_rarefation.
        """
        
        s = None
        
        if data_type in ['v', 'd', 'j', 'cdr3_prot_len', 'cdr3_nucl_len', 'v_del_len', 'j_del_len', 'vj_ins_len']:
            s = self.__get_freq(data_type, prob)
        elif data_type == 'vdj':
            s = self.__get_freq_vdj(data_type, prob)
        elif data_type == 'vdj_rarefaction':
            s = self.__get_est_vdj_rarefaction(data_type, func)
        else:
            raise ValueError(data_type + ' is not supported by get_summary method.')
        
        return s
    
    
    def __get_freq(self, data_type, prob = False):
        freq = None
        if data_type == 'v_del_len':
            freq = self.indels.v_del_len.value_counts(dropna = False)
            freq.index = ['ambiguous' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
            freq = freq.sort_index(ascending = True)
        elif data_type == 'j_del_len':
            freq = self.indels.j_del_len.value_counts(dropna = False)
            freq.index = ['ambiguous' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
            freq = freq.sort_index(ascending = True)
        elif data_type == 'vj_ins_len':
            freq = self.indels.vj_ins_len.value_counts(dropna = False)
            freq.index = ['ambiguous' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
            freq = freq.sort_index(ascending = True)
        elif data_type == 'cdr3_prot_len':
            freq = self.cdr3.prot_len.value_counts(dropna = False)
            freq.index = ['ambiguous' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
            freq = freq.sort_index(ascending = True)
        elif data_type == 'cdr3_nucl_len':
            freq = self.cdr3.nucl_len.value_counts(dropna = False)
            freq.index = ['ambiguous' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
            freq = freq.sort_index(ascending = True)
        elif data_type == 'v':
            freq = self.vdj.v.value_counts(dropna = False)
        elif data_type == 'd':
            freq = self.vdj.d.value_counts(dropna = False)
            freq.index = ['ambiguous' if (type(_i) == np.float and np.isnan(_i)) else _i for _i in freq.index]
        elif data_type == 'j':
            freq = self.vdj.j.value_counts(dropna = False)
        if prob:
            freq = freq / freq.sum()
        
        return freq
    
    
    
    def __get_freq_vdj(self, data_type, prob):
        __sep = '__________'
        vdj_combinations = self.vdj.v.replace(np.nan, 'NA') + __sep + \
                           self.vdj.d.replace(np.nan, 'NA') + __sep + \
                           self.vdj.j.replace(np.nan, 'NA')
        freq = vdj_combinations.value_counts(dropna = False)
        
        if prob:
            freq = freq / freq.sum(axis = 2)
        
        vdj_v = []
        vdj_d = []
        vdj_j = []
        for vdj_combination in freq.index:
            vv, dd, jj = vdj_combination.split(__sep)
            vdj_v.append(vv)
            vdj_d.append(dd)
            vdj_j.append(jj)
        freq = pd.DataFrame({'v': vdj_v, 'd': vdj_d, 'j': vdj_j, 'frequency': freq.values},
                             columns = ['v', 'd', 'j', 'frequency']).replace('NA', np.nan)
        return freq
    
    
    
    def __get_est_vdj_rarefaction(self, data_type, func = 'mean'):
        rdata = None
        if self.div.rarefaction['vdj'] is not None:
            if func == 'raw':
                rdata = self.div.rarefaction['vdj']
            elif func == 'mean':
                rdata = self.div.rarefaction['vdj'].mean(axis = 1)
        return rdata
 
        
        
    
    def samplingresampling_study(self, data = None, n = 1000):
        """Sampling-resampling study.
        
        Args:
            data (str): ``vdj`` or ``cdr3``. 
        
        Sampling-resamplig study.
        """
        
        if data is None or data == 'all':
            data = ['vdj']
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
        """Perform rarefaction study for VDJ combination or CDR3 sequenece.
        
        Args:
            data (str): One of ``vdj`` or ``cdr3`` can be specified for diversity study.
            n (int): The number of performing of capture-recapture procedures.
        """
        
        if data is None or data == 'all':
            data = ['vdj']
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
        """PyDAIRStatsRecords class initialize method.
        
        """
        
        self.__records = []
        self.__index = 0
    
    def __len__(self):
        return len(self.__records)
    
    def len(self):
        """Retrive the number of PyDAIRStatsRecord class object.
        
        Retrive the number of PyDAIRStatsRecord objects, thus the number of samples.
        """
        
        return self.__len__()

    def __iter__(self):
        self.__index = 0
        return self
    
    def __next__(self):
        try:
            stats_record = self.__records[self.__index]
        except IndexError:
            raise StopIteration
        self.__index += 1
        return stats_record
    
    def next(self):
        """Return the next IgSeq object from the iterator of PyDAIRStatsRecords class objects.
        
        """
        
        return self.__next__()
    
    
    def append(self, stats_record):
        """Append PyDAIRStatsRecord into PyDAIRStatsRecords class object.
        
        """
        
        self.__records.append(stats_record)
    
    
    def get_record(self, i):
        """Retrive A PyDAIRStatsRecord.
        
        Args:
            i (int): The index of PyDAIRStatsRecord should be returned.
        
        Returns:
            PyDAIRStatsRecord: Return a PyDAIRStatsRecord object.
        """
        
        try:
           stats_record = self.__records[i]
        except IndexError:
            raise IndexError
        return stats_record
    
    
    def set_record(self, i, stats_record):
        """Update PyDAIRStatsRecord into PyDAIRStatsRecords class object.
        
        Args: 
            i (int): The index of PyDAIRStatsRecord should be updated.
            stats_record: A PyDAIRStatsRecord that will used for updating.
        """
        
        self.__records[i] = stats_record


    
class PyDAIRStats:
    '''
    Statistics analysis of repertoire sequences.
    
    The class for storing the PyDAIR data file. Before use initialize this class,
    one should run create PyDAIR format files by other PyDAIR functions. This class
    stored many data of analyzed data.
    '''
    
    def __init__(self, pydair_file, pydair_id = None,
                 discard_ambiguous_D = False, productive_only = False):
        """PyDAIRStats class initialize method.
        
        Args:
            pydair_file (list): A list the contains multiple PYDAIR file path. 
            pydair_id (list): A list of sample names.
            discard_ambiguous_D (bool): If true, discard ambiguous D before analysis.
            productive_only (bool): If true, analyze sequences with stop codons.
        """
        
        if pydair_id is None:
            pydair_id = []
            for i in range(len(pydair_file)):
                pydair_id.append('individual ' + str(i + 1))
        
        self.__pydair_file   = pydair_file
        self.__pydair_id     = pydair_id
        self.__discard_ambiguous_D = discard_ambiguous_D
        self.__productive_only = productive_only
        self.samples   = None
        
        # parse PyDAIR files
        self.__parse_pydair_files()
        
    
    def __parse_pydair_files(self):
        self.samples  = PyDAIRStatsRecords()
        for i in range(len(self.__pydair_file)):
            v = []
            d = []
            j = []
            orf = []
            cdr3_prot_seq = []
            cdr3_nucl_seq = []
            v_del = []
            j_del = []
            vj_ins = []
            
            pydair_fh = PyDAIRIO(self.__pydair_file[i], 'r')
            for igseq in pydair_fh.parse():
                v.append(igseq.v.sbjct.name)
                d.append(igseq.d.sbjct.name)
                j.append(igseq.j.sbjct.name)
                orf.append(igseq.query.orf)
                
                cdr3_data = igseq.get_cdr3_data()
                if cdr3_data.nucl_seq is None:
                    cdr3_data.nucl_seq = ''
                    cdr3_data.prot_seq = ''
                if '*' in cdr3_data.prot_seq:
                    cdr3_data.nucl_seq = ''
                    cdr3_data.prot_seq = ''
                cdr3_prot_seq.append(cdr3_data.prot_seq)
                cdr3_nucl_seq.append(cdr3_data.nucl_seq)
                
                v_del.append(igseq.indels.v_deletion)
                j_del.append(igseq.indels.j_deletion)
                vj_ins.append(igseq.indels.vj_insertion)
                
            sample_record = PyDAIRStatsRecord(self.__pydair_id[i], v, d, j, orf,
                                              cdr3_nucl_seq, cdr3_prot_seq,
                                              v_del, j_del, vj_ins,
                                              self.__discard_ambiguous_D, self.__productive_only)
            self.samples.append(sample_record)
            pydair_fh.close()
    
    
    def rarefaction_study(self, data = None, n = 1000):
        """Rarefaction study for estimating the number of VDJ combination.
        
        Args:
            data (str): 'vdj' or 'cdr3'.
        
        Rarefaction study for estimating the number of VDJ combination for each sample.
        """
        
        for bsample in self.samples:
            bsample.rarefaction_study(data, n)
    
    
    
    def get_summary(self, data_type, prob = False):
        '''Get the frequence of usage with data frame class object.
        
        Args:
            data_type (str): A string to specify data type. ``v``, 
                             ``d``, ``j``, ``vdj``, ``cdr3_nucl_len``,
                             ``cdr3_prot_len``, ``v_del_len``, ``j_del_len``,
                             and ``vj_ins_len`` are supported.
            freq (bool): If `True`, return the frequences of counts.
            prob (bool): If prob is `True`, return the probability, and omit `freq`.
        
        Returns:
            A Pandas DataFrame class object.
        
        Get the frequence of usages of the ``data_type`` as DataFrame class object.
        '''
        
        s = None
        
        if data_type in ['v', 'd', 'j', 'cdr3_prot_len', 'cdr3_nucl_len', 'v_del_len', 'j_del_len', 'vj_ins_len']:
            s = self.__get_freq(data_type, prob)
        elif data_type == 'vdj':
            s = self.__get_freq_vdj(data_type, prob)
        elif data_type == 'vdj_rarefaction':
            s = self.__get_est_vdj_rarefaction(data_type)
        else:
            raise ValueError(data_type + ' is not supported by get_summary method.')
        
        return s
        
        
        
    def __get_freq(self, data_type, prob = False):
        sample_freqs = []
        sample_names = []
        for sample in self.samples:
            sample_freqs.append(sample.get_summary(data_type, prob = prob))
            sample_names.append(sample.name)
        freq_dataframe = pd.concat(sample_freqs, axis = 1)
        freq_dataframe.columns = sample_names
        
        if prob:
            freq_dataframe = freq_dataframe / freq_dataframe.sum(axis = 1)
        
        freq_dataframe.columns = sample_names
        
        if data_type in ['v', 'd', 'j']:
            freq_dataframe = freq_dataframe.ix[freq_dataframe.mean(axis = 1).sort_values(ascending = False).index]
        else:
            freq_dataframe = freq_dataframe.set_index([[int(dx) for dx in freq_dataframe.index.values]])
        return freq_dataframe
    
    
    def __get_freq_vdj(self, data_type, prob = False):
        __sep = '__________'
        sample_freqs = []
        sample_names = []
        for sample in self.samples:
            freq = sample.get_summary(data_type, prob = prob)
            freqval = pd.Series(freq.frequency)
            freqval.index = freq.v.replace(np.nan, 'NA') + __sep + \
                            freq.d.replace(np.nan, 'NA') + __sep + \
                            freq.j.replace(np.nan, 'NA')
            sample_freqs.append(freqval)
            sample_names.append(sample.name)
        freq_dataframe = pd.concat(sample_freqs, axis = 1)
        freq_dataframe.columns = sample_names
        
        vdj_v = []
        vdj_d = []
        vdj_j = []
        for vdj_combination in freq_dataframe.index:
            vv, dd, jj = vdj_combination.split(__sep)
            vdj_v.append(vv)
            vdj_d.append(dd)
            vdj_j.append(jj)
        
        freq_vdjcmbn = pd.concat([pd.Series(vdj_v), pd.Series(vdj_d), pd.Series(vdj_j)], axis = 1).replace('NA', np.nan)
        freq_vdjcmbn.columns = ['V', 'D', 'J']
        
        freq_vdjcmbn.reset_index(drop = True, inplace = True)
        freq_dataframe.reset_index(drop = True, inplace = True)
        
        freq = pd.concat([freq_vdjcmbn, freq_dataframe], axis = 1)
        return freq
    
    
    
    def __get_est_vdj_rarefaction(self, data_type):
        rdata = []
        sample_names = []
        for sample in self.samples:
            sample_rdata = sample.get_summary('vdj_rarefaction', 'mean')
            sample_names.append(sample.name)
            if sample_rdata is not None:
                rdata.append(sample_rdata)
        if len(rdata) > 0:
            rdata = pd.concat(rdata, axis = 1)
            rdata.columns = sample_names
        else:
            rdata = None
        
        return rdata
        
    
    
