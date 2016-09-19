from __future__ import division
import os
import sys
import warnings
import random
import re
import argparse
import subprocess
import logging
from PyDAIR import *
from PyDAIR.io.PyDAIRIO import *
from PyDAIR.utils.PyDAIRArgs import *
from PyDAIR.app.PyDAIRAPP import *
from PyDAIR.stats.PyDAIRStats import *
from PyDAIR.plot.PyDAIRPlot import *

from Bio import SeqIO


logging.basicConfig(level = logging.INFO, format = '%(levelname)-8s %(message)s')


class PyDAIRAPPParseSeq:
    '''
    Framework of PyDAIR software
    
    This class provides the first phase analysis. The parameters should be given by
    PyDAIRArgs class object.
    '''
    
    def __init__(self, args):
        '''Initialize PyDAIR class object.
        
        Fill up all parameters for esecuting PyDAIR tool at first.
        '''
        self.__args = args
        # file names used in PyDAIR
        self.__pydair_output_vj  = self.__args.output_path + '.vj.pydair'
        self.__pydair_output_vdj = self.__args.output_path + '.vdj.pydair'
        
        self.__blast_output_prefix = self.__args.output_path
        self.__blast_output_v = self.__blast_output_prefix + '.v.blast.tsv'
        self.__blast_output_d = self.__blast_output_prefix + '.d.blast.tsv'
        self.__blast_output_j = self.__blast_output_prefix + '.j.blast.tsv'
        self.__unaligned_fasta = self.__blast_output_prefix + '.unaligned.fa'
        self.__cdr3_fasta      = self.__blast_output_prefix + '.cdr3.fa'
        
        # data
        self.__pydair_records = []
        
        self.__species = self.__args.species
        
        # protocol log
        self.__log_parsed_v = False
        self.__log_parsed_d = False
        self.__log_parsed_j = False
        self.__log_blasted_v = False
        self.__log_blasted_d = False
        self.__log_blasted_j = False
        
        
    
    def blast(self, gene = None, program_path = ''):
        '''Run BLAST.
        
        Use 'blast' command to assign V, D, or J genes.
        Note that the 'blast' command should be exported.
        '''
        # parameter settings
        match_score, mismatch_score, gap_open_penalty, gap_extend_penalty, \
        word_size, evalue_cutoff, blast_db, \
        query_file, blast_output_path = self.__get_blast_params(gene)
        
        if os.path.exists(blast_output_path):
            logging.info('Found BLAST results for ' + gene + ' genes (' + blast_output_path + '), use it instead of executing BLAST.')
        else:
            # BLAST command
            blast_cmd = '''
blastn -db %s -query %s -out %s -word_size %s -reward %s -penalty %s -gapopen %s -gapextend %s -outfmt "6 qseqid qstart qend qseq sseqid sstart send sseq score pident sstrand evalue" -num_alignments 1
'''
            blast_cmd = blast_cmd % (blast_db, query_file, blast_output_path,
                                     word_size, match_score, mismatch_score,
                                     gap_open_penalty, gap_extend_penalty)
            blast_cmd = blast_cmd.rstrip()
            # perform BLAST
            logging.info('PyDAIR: Start to perform BLAST for ' + gene + ' gene.' )
            logging.info(blast_cmd)
            proc_sig = subprocess.call(program_path + ' ' + blast_cmd, shell = True)
            logging.info('PyDAIR: Finished performing BLAST for ' + gene + ' gene.' )
        
        if gene == 'v':
            self.__log_blasted_v = True
        elif gene == 'd':
            self.__log_blasted_d = True
        elif gene == 'j':
            self.__log_blasted_j = True
    
        
    def __create_dict_from_fasta(self, fasta_path):
        fa_dict = {}
        with open(fasta_path, 'r') as fa_fh:
            for record in SeqIO.parse(fa_fh, 'fasta'):
                fa_dict[record.id] = str(record.seq).lower()
        return fa_dict
    

    def __get_blast_params(self, gene = None):
        match_score        = None
        mismatch_score     = None
        gap_open_penalty   = None
        gap_extend_penalty = None
        word_size          = None
        evalue_cutoff      = None
        blast_db           = None
        query_file         = None
        blast_output_path  = None
        
        if gene == 'v':
            match_score        = self.__args.v_align_args.match
            mismatch_score     = self.__args.v_align_args.mismatch
            gap_open_penalty   = self.__args.v_align_args.gapopen
            gap_extend_penalty = self.__args.v_align_args.gapextend
            word_size          = self.__args.v_align_args.wordsize
            evalue_cutoff      = self.__args.v_align_args.cutoff
            blast_db           = self.__args.v_align_args.db
            query_file         = self.__args.q_file_path
            blast_output_path  = self.__blast_output_v
        elif gene == 'd': 
            match_score        = self.__args.d_align_args.match
            mismatch_score     = self.__args.d_align_args.mismatch
            gap_open_penalty   = self.__args.d_align_args.gapopen
            gap_extend_penalty = self.__args.d_align_args.gapextend
            word_size          = self.__args.d_align_args.wordsize
            evalue_cutoff      = self.__args.d_align_args.cutoff
            blast_db           = self.__args.d_align_args.db
            query_file         = self.__unaligned_fasta
            blast_output_path  = self.__blast_output_d
        elif gene == 'j':
            match_score        = self.__args.j_align_args.match
            mismatch_score     = self.__args.j_align_args.mismatch
            gap_open_penalty   = self.__args.j_align_args.gapopen
            gap_extend_penalty = self.__args.j_align_args.gapextend
            word_size          = self.__args.j_align_args.wordsize
            evalue_cutoff      = self.__args.j_align_args.cutoff
            blast_db           = self.__args.j_align_args.db
            query_file         = self.__args.q_file_path
            blast_output_path  = self.__blast_output_j
        else:
            raise ValueError('The argument should be one of "v", "d" and "j".')
        return [match_score, mismatch_score, gap_open_penalty, gap_extend_penalty,
                word_size, evalue_cutoff, blast_db, query_file, blast_output_path]
        
    
    def parse_VJ(self):
        '''Parse BLAST results of V and J.
        
        This method should run after having the BLAST results of V and J.
        This method performs:
            1. parse BLAST result of V gene, identify the used V gene.
            2. parse BLAST result of V gene, identify the used J gene.
            3. save the parsed results into self.__pydair_records.
        '''
        # load FASTA file into dictionary
        fa_q_dict = self.__create_dict_from_fasta(self.__args.q_file_path)
        fa_v_dict = self.__create_dict_from_fasta(self.__args.v_file_path)
        fa_j_dict = self.__create_dict_from_fasta(self.__args.j_file_path)
        # create dictionary from BLAST results
        vj_dict =  {}
        self.__create_dict_from_blastoutput(self.__blast_output_v, 'v', vj_dict, self.__args.v_align_args.cutoff)
        self.__create_dict_from_blastoutput(self.__blast_output_j, 'j', vj_dict, self.__args.j_align_args.cutoff)
        
        # parse BLAST results
        for seq_name, gene_dict in vj_dict.items():
            ## 0:qseqid 1:qstart 2:qend 3:qseq 4:sseqid 5:sstart 6:send 7:sseq 8:score 9:pident 10:sstrand"
            v = gene_dict['v']
            j = gene_dict['j']
            igseqv = None
            igseqd = None
            igseqj = None
            # if not None, set values
            if v != None:
                igseqv_q = IgSeqAlignQuery(seq_name, fa_q_dict[seq_name], v[3], v[1], v[2], '+')
                igseqv_s = IgSeqAlignSbjct(    v[4], fa_v_dict[v[4]], v[7], v[5], v[6], v[10])
                igseqv = IgSeqAlign(igseqv_q, igseqv_s, v[9], v[11])
            if j != None:
                if (v == None) or (v != None and v[10] == j[10]):
                    igseqj_q = IgSeqAlignQuery(seq_name, fa_q_dict[seq_name], j[3], j[1], j[2], '+')
                    igseqj_s = IgSeqAlignSbjct(    j[4], fa_j_dict[j[4]], j[7], j[5], j[6], j[10])
                    igseqj = IgSeqAlign(igseqj_q, igseqj_s, j[9], j[11])
            # create IgSeq class object
            igseq = IgSeq(self.__species, igseqv, igseqd, igseqj)
            # find CDR3 region
            igseq.seek_cdr3()
            # if alignemnt of V and J are correct, then print out the results into file
            if igseq.valid_alignment:
                self.__pydair_records.append(igseq)
        
        self.__log_parsed_v = True
        self.__log_parsed_j = True
    
        
    
    def __create_dict_from_blastoutput(self,  blast_output_path, sbjct_db_name, parsed_data_dict, cutoff):
        with open(blast_output_path, 'r') as fh:
            for record in fh:
                record = record.replace('\n', '')
                cols = record.split('\t')
                if float(cols[11]) < cutoff:
                    # create new key if not defined
                    if cols[0] not in parsed_data_dict:
                        parsed_data_dict[cols[0]] = {'v': None, 'j': None, 'd': None}
                    # change plus/minus to +/-
                    if cols[10] == 'plus':
                        cols[10] = '+'
                    elif cols[10] == 'minus':
                        cols[10] = '-'
                    if parsed_data_dict[cols[0]][sbjct_db_name] == None:
                        parsed_data_dict[cols[0]][sbjct_db_name] = cols

    
    
    def write_pydair(self, file_name = None, file_format = 'pydair'):
        
        if file_name is None:
            if self.__log_parsed_v and self.__log_parsed_j and self.__log_parsed_d:
                file_name = self.__pydair_output_vdj
            elif self.__log_parsed_v and self.__log_parsed_j:
                file_name = self.__pydair_output_vj
        else:
            if self.__log_parsed_v and self.__log_parsed_j and self.__log_parsed_d:
                self.__pydair_output_vdj = file_name
            elif self.__log_parsed_v and self.__log_parsed_j:
                self.__pydair_output_vj = file_name
        
        if file_format == 'simple':
            file_name = file_name + '.simple'
        
        # open file and write Igseq object into file
        pydair_o_fh = PyDAIRIO(file_name, 'w', file_format)
        # read record in PyDAIR flat file
        for igseq in self.__pydair_records:
            pydair_o_fh.write(igseq)
        pydair_o_fh.close()
        logging.info('Data has been saved into ' + file_name + '.\n')
    
    
    
    def write_fasta(self, seq_type = 'unaligned_seq', mol_type = 'nucl'):
        '''Create unaligned seuqence FASTA
        
        '''
        if seq_type == 'unaligned_seq':
            self.__write_fasta_unaligned_seq(self.__unaligned_fasta)
        elif seq_type == 'cdr3':
            self.__write_fasta_cdr3_seq(self.__cdr3_fasta, mol_type)
        else:
           raise ValueError('seq_type should be unaligned.')
    
    
    def __write_fasta_cdr3_seq(self, filename, mol_type):
        with open(filename, 'w') as output_fa_fh:
            if os.path.exists(self.__pydair_output_vdj):
                pydair_fh = PyDAIRIO(self.__pydair_output_vdj, 'r', self.__args.pydair_format)
            else:
                pydair_fh = PyDAIRIO(self.__pydair_output_vj, 'r', self.__args.pydair_format)
            
            for igseq in pydair_fh.parse():
                cdr3_data = igseq.get_cdr3_data()
                if mol_type == 'nucl' and (cdr3_data.nucl_seq is not None and cdr3_data.nucl_seq != ''):
                    output_fa_fh.write('>' + igseq.query.name + '\n' + \
                                       cdr3_data.nucl_seq + '\n')
                elif mol_type == 'prot' and (cdr3_data.prot_seq is not None and cdr3_data.prot_seq != ''):
                    output_fa_fh.write('>' + igseq.query.name + '\n' + \
                                       cdr3_data.prot_seq + '\n')
        logging.info('Data (cdr3 sequences) has been saved into ' + filename + '.\n')
    
    
    def __write_fasta_unaligned_seq(self, filename):
        '''Write unaligned sequence FASTA
        
        If CDR3 is identified, use CDR3 region as unaligned sequence,
        otherwise use original un-aligned sequence.
        
        '''
        with open(filename, 'w') as output_fa_fh:
            if os.path.exists(self.__pydair_output_vdj):
                pydair_fh = PyDAIRIO(self.__pydair_output_vdj, 'r', self.__args.pydair_format)
            else:
                pydair_fh = PyDAIRIO(self.__pydair_output_vj, 'r', self.__args.pydair_format)
            
            for igseq in pydair_fh.parse():
                # first check CDR3, if CDR3 is not exsisted, then use unaligned region
                if igseq.variable_region.cdr3 is not None:
                    output_fa_fh.write('>' + igseq.query.name  + ' [CDR3_region]\n' + \
                        igseq.query.seq[(igseq.variable_region.cdr3[0] - 1):(igseq.variable_region.cdr3[1] - 1)] + '\n')
                else:
                    if igseq.variable_region.untemplate_region is not None:
                        untemplate_region_seq = igseq.query.seq[(igseq.variable_region.untemplate_region[0] - 1):(igseq.variable_region.untemplate_region[1] - 1)]
                        if untemplate_region_seq != '':
                            output_fa_fh.write('>' + igseq.query.name  + ' [Unaligned_region]\n' + untemplate_region_seq + '\n')
        logging.info('Data (unaligned sequences) has been saved into ' + filename + '.\n')
    
    
    def parse_VDJ(self):
        '''Parse BLAST results of D.
        
        This method should run after 'parse_VDJ' and after having the BLAST results of D.
        This method performs:
            1. parse BLAST results of D gene, identify the used D gene.
            2. reading the results of PyDAIR file (VJ results), and add 
        '''
        # create dictionary from FASTA file (only D gene)
        fa_d_dict = self.__create_dict_from_fasta(self.__args.d_file_path)
        # create dictionary from BLAST reuslts
        d_dict = {}
        self.__create_dict_from_blastoutput(self.__blast_output_d, 'd', d_dict, self.__args.d_align_args.cutoff)
        # save to PyDAIR format file with VJ results
        pydair_vj  = PyDAIRIO(self.__pydair_output_vj,  'r', self.__args.pydair_format)
        for i in range(len(self.__pydair_records)):
            igseq = self.__pydair_records[i]
            igseqd = None
            if igseq.query.name in d_dict:
                d = d_dict[igseq.query.name]['d']
                # the strands of V and J genes are same at this step
                if d != None and igseq.v.sbjct.strand == d[10]:
                    igseqd_q = IgSeqAlignQuery(igseq.query.name, igseq.query.seq, d[3], d[1], d[2], '+')
                    igseqd_s = IgSeqAlignSbjct(            d[4], fa_d_dict[d[4]], d[7], d[5], d[6], d[10])
                    igseqd = IgSeqAlign(igseqd_q, igseqd_s, d[9], d[11])
                    igseq.set_igseqalign(igseqd, 'd')
            self.__pydair_records[i] = igseq
        self.__log_parsed_d = True









class PyDAIRAPPStats:
    '''
    
    '''
    def __init__(self, args):
        self.__args = args
        self.__sample_names        = args.sample_names
        self.__pydair_files         = args.pydair_files
        self.__contain_ambiguous_D = args.contain_ambiguous_D
        self.__contain_stopcodon   = args.contain_stopcodon
        self.__output_prefix       = args.output_prefix
        self.__figure_format       = args.figure_format
        self.__out_file_rarefaction  = self.__output_prefix + '.rarefaction.tsv'
        self.__out_file_samplingresamplig  = self.__output_prefix + '.samplingresamplig.tsv'
        
        # Data files
        self.__o_file_v_freq = []
        self.__o_file_d_freq = []
        self.__o_file_j_freq = []
        self.__o_file_vdj_freq = []
        self.__o_file_rarefaction = []
        self.__o_file_samplingresampling = []
        self.__o_file_cdr3_prot_len_freq = []
        self.__o_file_cdr3_nucl_len_freq = []
        for sample_name in self.__sample_names:
            sample_name_f = sample_name.replace(' ', '_')
            self.__o_file_v_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.v.freq.tsv')
            self.__o_file_d_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.d.freq.tsv')
            self.__o_file_j_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.j.freq.tsv')
            self.__o_file_vdj_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.vdj.freq.tsv')
            self.__o_file_rarefaction.append(self.__output_prefix + '.sample.' + sample_name_f + '.rarefaction.tsv')
            self.__o_file_samplingresampling.append(self.__output_prefix + '.sample.' + sample_name_f + '.samplingresampling.tsv')
            self.__o_file_cdr3_prot_len_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.cdr3_prot_length.freq.tsv')
            self.__o_file_cdr3_nucl_len_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.cdr3_nucl_length.freq.tsv')
        
        # Plot files
        self.__o_figure_v_freq   = self.__output_prefix + '.v.freq.' + self.__figure_format
        self.__o_figure_d_freq   = self.__output_prefix + '.d.freq.' + self.__figure_format
        self.__o_figure_j_freq   = self.__output_prefix + '.j.freq.' + self.__figure_format
        self.__o_figure_vdj_freq = []
        for sample_name in self.__sample_names:
            sample_name_f = sample_name.replace(' ', '_')
            self.__o_figure_vdj_freq.append(self.__output_prefix + '.sample.' + sample_name_f + '.vdj.freq.' + self.__figure_format)
        
        self.__o_figure_rarefaction        = self.__output_prefix + '.vdj.rarefaction.' + self.__figure_format
        self.__o_figure_samplingresampling = self.__output_prefix + '.cdr3.samplingresampling.' + self.__figure_format
        self.__o_figure_cdr3_len_dist      = self.__output_prefix + '.cdr3.length.dist.' + self.__figure_format
        
        self.stats = PyDAIRStats(self.__pydair_files, 'pydair', self.__sample_names,
                                 self.__contain_ambiguous_D, self.__contain_stopcodon)
        self.plots = PyDAIRPlot(self.stats)
    
    
    
    
    def write_cdr3_len_freq(self):
        '''
        Write the length distribution of CDR3 sequence.
        
        Write the length distribution of CDR3 nucleotide and protein sequences. The output file
        is tab-separator text file, and consists of two columns. The first column is the length
        of CDR3 sequence, and the second column is frequency related the length.
        '''
        _freq_nucl = []
        _freq_prot = []
        for sample_i in range(len(self.stats.samples)):
            bsample = self.stats.samples.get_record(sample_i)
            _freq_nucl.append(bsample.get_freq('cdr3_nucl_len'))
            _freq_prot.append(bsample.get_freq('cdr3_prot_len'))
        
        for _i in range(len(self.stats.samples)):
            _freq_nucl[_i].to_csv(self.__o_file_cdr3_nucl_len_freq[_i], sep = '\t', na_rep = 'NA',
                                  header = True, index = True)
            _freq_prot[_i].to_csv(self.__o_file_cdr3_prot_len_freq[_i], sep = '\t', na_rep = 'NA',
                                  header = True, index = True)
    
    
    
    
    def write_diversity(self, methods = ['rarefaction', 'samplingresampling']):
        '''
        Write the results of diversity studies into text file.
        
        Write the reuslts of diversity studies into tab-separator files. The file has multiple
        columns. The first column indicates that the smapling size, and from the second columns
        indicate the results of sampling-resampling study.
        '''
        if not isinstance(methods, list):
            methods = [methods]
        
        for method in methods:
            for sample_i in range(len(self.stats.samples)):
                bsample = self.stats.samples.get_record(sample_i)
                x = None
                y = None
                file_path = None
                
                if method == 'rarefaction':
                    if bsample.rarefaction is not None:
                        x = bsample.rarefaction.x
                        y = bsample.rarefaction.y
                        file_path = self.__o_file_rarefaction[sample_i]
                elif method == 'samplingresampling':
                    if bsample.samplingresampling is not None:
                        x = bsample.rarefaction.x
                        y = bsample.rarefaction.y
                        file_path = self.__o_file_samplingresampling[sample_i]
                else:
                    raise ValueError('The \'methods\' argument of \'write_diversity\' method should be one of \'rarefaction\' and \'samplingresampling\'.')
                
                if x is None and y is None:
                    continue
                dat = []
                for x_i in range(len(x)):
                    _dat = [x[x_i]]
                    for y_i in range(len(y)):
                        _dat.append(y[y_i][x_i])
                    dat.append(_dat)
                with open(file_path, 'w') as fh:
                    for i in range(len(dat)):
                        fh.write('\t'.join([str(_dat).replace('None', 'NA') for _dat in dat[i]]) + '\n')
    
    
    
    
    
    def write_freq(self,  genes = ['v', 'd', 'j', 'vdj']):
        '''
        Write frequency of V, D, J genes and VDJ combinations into text file.
        The output file is tab-separator format text file. If 'genes' option is one of
        'v', 'd', 'j', the output file consists of two columns, and the first column
        is gene name and the second column is the frequency of the related gene. If
        'genes' option is 'vdj', the output file consists of four columns, the first
        three columns is the gene names of V, D, and J, and the fifth column is the
        frequency of the related combination of V, D and J.
        '''
        if not isinstance(genes, list):
            genes = [genes]
        
        for gene in genes:
            for sample_i in range(len(self.stats.samples)):
                bsample = self.stats.samples.get_record(sample_i)
                _freq = bsample.get_freq(gene)
                
                if gene == 'v' or gene == 'V':
                    file_path = self.__o_file_v_freq[sample_i]
                elif gene == 'd' or gene == 'D':
                    file_path = self.__o_file_d_freq[sample_i]
                elif gene == 'j' or gene == 'J':
                    file_path = self.__o_file_j_freq[sample_i]
                elif gene == 'vdj' or gene == 'VDJ':
                    file_path = self.__o_file_vdj_freq[sample_i]
                else:
                    raise ValueError('The \'genes\' argument of \'write_freq\' method should be one of \'v\', \'d\', \'j\', and \'vdj\'.')
                _freq.to_csv(file_path, sep = '\t', na_rep = 'NA', header = True, index = True)
    
    
    
    
    def plot(self):
        self.plots.plot_freq(gene = 'v', fig_name = self.__o_figure_v_freq)
        self.plots.plot_freq(gene = 'd', fig_name = self.__o_figure_d_freq)
        self.plots.plot_freq(gene = 'j', fig_name = self.__o_figure_j_freq)
        self.plots.plot_freq(gene = 'vdj', fig_name = self.__o_figure_vdj_freq)
        self.plots.plot_rarefaction(fig_name = self.__o_figure_rarefaction)
        self.plots.plot_samplingresampling(fig_name = self.__o_figure_samplingresampling)
        self.plots.plot_cdr3_len_dist(fig_name = self.__o_figure_cdr3_len_dist)
    
    
