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
#from PyDAIR.plot.PyDAIRPlot import *

from Bio import SeqIO


logging.basicConfig(level = logging.INFO, format = '%(levelname)-8s %(message)s')


class PyDAIRAPPParseSeq:
    """PyDAIRAPPParseSeq class.
    
    This class implements the methods to parse Rep-Seq sequencces.
    It provides functions of
    (i) V, D, and J segment identification;
    (ii) CDR3 segment identificaiton;
    (iii) 3'-end V deletion detection;
    (iv) 5'-end J deletion detection;
    and (v) insertion during VDJ junctions detection.
    """
    
    def __init__(self, args):
        """PyDAIRAPPParseSeq class initialize method.
        
        Args:
            args (PyDAIRParseSeqArgs): PyDAIRParseSeqArgs class object that contains 
                                       path to Rep-Seq data and BLAST parameters.
        
        >>> valign = PyDAIRBlastArgs('balstdb_v', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 21, eval_cutoff = 1e-90)
        >>> dalign = PyDAIRBlastArgs('balstdb_d', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 2, eval_cutoff = 1e-3)
        >>> jalign = PyDAIRBlastArgs('balstdb_j', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 9, eval_cutoff = 1e-10)
        >>> args = PyDAIRParseSeqArgs('fugu', q = 'q.fa', v = 'v.fa', d = 'd.fa',
        >>>                           j = 'j.fa', o = 'out_prefix', f = 'pydair',
        >>>                           valign = valign, dalign = dalign, jalign = jalign)
        >>> obj = PyDAIRAPPParseSeq(args)
        """
        
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
        
        # motif
        self.__v_motif = self.__args.v_motif
        self.__j_motif = self.__args.j_motif
        
        # protocol log
        self.__log_parsed_v = False
        self.__log_parsed_d = False
        self.__log_parsed_j = False
        self.__log_blasted_v = False
        self.__log_blasted_d = False
        self.__log_blasted_j = False
        
        
    
    def blast(self, gene = None, program_path = ''):
        """Execute BLAST.
        
        This method internally calls `blastn` program to align Rep-Seq sequences
        against the germeline (**gene**) database.
        The BLAST result is saved in the specified output directory with TSV format.
        The BLAST parameters are given by PyDAIRBlastArgs class object through
        `__get_blast_params` method.
        
        Args:
            gene (str): Specified `v`, `d`, or `j` for executing BLAST.
            program_path (str): A path to `blastn` programa.
        
        """
        
        # parameter settings
        match_score, mismatch_score, gap_open_penalty, gap_extend_penalty, \
        word_size, evalue_cutoff, blast_db, \
        query_file, blast_output_path = self.__get_blast_params(gene)
        
        if os.path.exists(blast_output_path):
            logging.info('Found BLAST results for ' + gene + ' genes (' + blast_output_path + '), use it instead of executing BLAST.')
        else:
            # BLAST command
            blast_cmd = """
blastn -db %s -query %s -out %s -word_size %s -reward %s -penalty %s -gapopen %s -gapextend %s -outfmt "6 qseqid qstart qend qseq sseqid sstart send sseq score pident sstrand evalue" -num_alignments 1
"""
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
        """Create a dictionary class object from FASTA file.
        
        Args:
            fasta_path (str):
                A path to FASTA file.
        
        Returns:
            A dictionary with the keys of sequence names and the values of sequences.
        
        Read FASTA file and store all data into a dictionary class object.
        The keys of the dictionary are the sequence names and the values are
        the nucelotide sequences.
        """
        
        fa_dict = {}
        with open(fasta_path, 'r') as fa_fh:
            for record in SeqIO.parse(fa_fh, 'fasta'):
                fa_dict[record.id] = str(record.seq).lower()
        return fa_dict
    

    def __get_blast_params(self, gene = None):
        """Retrieve the BLAST parameters for a germline gene.
        
        Args:
            gene (str):
                A gene name. One of 'v', 'd', and 'j' should be specified.
        
        Returns:
            A list class object contains BLAST paramters such as `match_score`, `mismatch_score`, and so on.
        
        This method is called by `blast` method, and return a list class object that
        contains BLAST parameters. The parameters should be set in `self` before
        use this method.
        """
        
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
        """Identify V and J segments.
        
        Identify V and J segments that used in Rep-Seq data according to the BLAST
        outputs of V and J genes.
        This method is performed as the three steps:
        (i) execute BLAST against V gene database and identify V segment;
        (ii) execute BLAST against J gene database and identify J segment;
        and (iii) save the results into **IgSeq** class object.
        
        """
        
        const_tag_obj = IgConstantTag(self.__v_motif, self.__j_motif)
        
        
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
            igseq = IgSeq(igseqv, igseqd, igseqj)
            # find CDR3 region
            igseq.seek_cdr3(const_tag_obj)
            igseq.seek_orf()
            igseq.find_indels()
            # if alignemnt of V and J are correct, then print out the results into file
            if igseq.valid_alignment:
                self.__pydair_records.append(igseq)
        
        self.__log_parsed_v = True
        self.__log_parsed_j = True
    
        
    
    def __create_dict_from_blastoutput(self,  blast_output_path, sbjct_db_name, parsed_data_dict, cutoff):
        """Identify one of V, D, and J genes and save reuslts into a dictionary class object.
        
        Args:
            blast_output_path (str):
                A path to BLAST result.
            sbjct_db_name (str):
                A gene name. One of 'v', 'd', and 'j' should be specified.
            parsed_data_dict (str):
                A two-level dictionary class object. The keys of first \
                dictionary are sequence names that in Rep-Seq data. \
                The keys of second dicitonary are 'v', 'd' and 'j'. \
                The vlaues of second dictionary are the assigned gene name \
                in BLAST database (i.e., identification results).
            cutoff (float):
                A threshold of E-value to cutoff BLAST hits.
        
        Returns:
            There is no return object. However, the `parsed_data_dict`
            is given by the parental methods and shares with parental methods.
        
        Read the BLAST results and find the best hit of each Rep-Seq reads,
        to identify V, D, or J genes.
        """
        
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

    
    
    def write_pydair(self, file_name = None):
        """Write results into file with PYDAIR format.
        
        Write parsed results into file with PYDAIR format.
        
        Args:
            file_name (str): A file path to save the result.
        
        """
        
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
        
        # open file and write Igseq object into file
        pydair_o_fh = PyDAIRIO(file_name, 'w')
        # read record in PyDAIR flat file
        for igseq in self.__pydair_records:
            pydair_o_fh.write(igseq)
        pydair_o_fh.close()
        logging.info('Data has been saved into ' + file_name + '.\n')
    
    
    
    def write_fasta(self, seq_type = 'unaligned_seq', mol_type = 'nucl'):
        """Write FASTA file.
        
        Write unaligned sequences or CDR3 sequences into specified file
        with FASTA format.
        
        Args:
            seq_type (str): A sequence type to save. One of `unaligned_seq`
                            or `cdr3` can be specified.
            mol_type (str): A molecular type. Only 'nucl' is supported.
        
        """
        
        if seq_type == 'unaligned_seq':
            self.__write_fasta_unaligned_seq(self.__unaligned_fasta)
        elif seq_type == 'cdr3':
            self.__write_fasta_cdr3_seq(self.__cdr3_fasta, mol_type)
        else:
           raise ValueError('seq_type should be unaligned.')
    
    
    def __write_fasta_cdr3_seq(self, filename, mol_type):
        """Write FASTA file of CDR3 sequences.
        
        Args:
            filename (str):
                A file path to save into FASTA file.
            mol_type (str):
                A file format to specified the molecular types in FASTA.
        
        Write CDR3 sequences with nucleotide (`nucl`) or amino acid (`prot`).
        """
        
        with open(filename, 'w') as output_fa_fh:
            if os.path.exists(self.__pydair_output_vdj):
                pydair_fh = PyDAIRIO(self.__pydair_output_vdj, 'r')
            else:
                pydair_fh = PyDAIRIO(self.__pydair_output_vj, 'r')
            
            for igseq in pydair_fh.parse():
                cdr3_data = igseq.get_cdr3_data()
                if mol_type == 'nucl' and (cdr3_data.nucl_seq is not None and cdr3_data.nucl_seq != ''):
                    output_fa_fh.write('>' + igseq.query.name + '\n' + \
                                       cdr3_data.nucl_seq + '\n')
                elif mol_type == 'prot' and (cdr3_data.prot_seq is not None and cdr3_data.prot_seq != ''):
                    output_fa_fh.write('>' + igseq.query.name + '\n' + \
                                       cdr3_data.prot_seq + '\n')
            pydair_fh.close()
        logging.info('Data (cdr3 sequences) has been saved into ' + filename + '.\n')
    
    
    def __write_fasta_unaligned_seq(self, filename):
        """Write unaligned sequence FASTA
        
        Args:
            filename (str):
                A path to save FASTA.
        
        If CDR3 is identified, use CDR3 region as unaligned sequence,
        otherwise use original un-aligned sequence.
        """
        
        with open(filename, 'w') as output_fa_fh:
            if os.path.exists(self.__pydair_output_vdj):
                pydair_fh = PyDAIRIO(self.__pydair_output_vdj, 'r')
            else:
                pydair_fh = PyDAIRIO(self.__pydair_output_vj, 'r')
            
            for igseq in pydair_fh.parse():
                if igseq.indels is not None and igseq.indels.vj_insertion is not None:
                    if len(igseq.indels.vj_insertion) > 4:
                        output_fa_fh.write('>' + igseq.query.name + '\n' + igseq.indels.vj_insertion + '\n')
            pydair_fh.close()
        logging.info('Data (unaligned sequences >= 5 nt.) has been saved into ' + filename + '.\n')
    
    
    def parse_VDJ(self):
        """Identify D segment.
        
        This method should run after **parse_VJ**.
        This method is performed as two steps:
        (i) execute BLAST against D gene database
        and (ii) identify D segment.
        """
        
        # create dictionary from FASTA file (only D gene)
        fa_d_dict = self.__create_dict_from_fasta(self.__args.d_file_path)
        # create dictionary from BLAST reuslts
        d_dict = {}
        self.__create_dict_from_blastoutput(self.__blast_output_d, 'd', d_dict, self.__args.d_align_args.cutoff)
        # save to PyDAIR format file with VJ results
        pydair_vj  = PyDAIRIO(self.__pydair_output_vj,  'r')
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
        pydair_vj.close()









class PyDAIRAPPStats:
    """PyDAIRAPPStats class.
    
    This class implements the methods to summarize parsed results.
    """
    
    def __init__(self, args):
        """PyDAIRAPPStats class initialize method.
        
        Args:
            args (PyDAIRStatsArgs): PyDAIRStatsArgs class object.
        """
        
        self.__args = args
        self.__sample_names        = args.sample_names
        self.__pydair_files        = args.pydair_files
        self.__discard_ambiguous_D = args.discard_ambiguous_D
        self.__productive_only   = args.productive_only
        self.__output_prefix       = args.output_prefix
        self.__estimate_vdj_combination = args.estimate_vdj_combination
        
        # data files
        __o_path = {
            'v': {},
            'd': {},
            'j': {},
            'vdj': {},
            'cdr3_nucl_len': {},
            'cdr3_prot_len': {},
            'v_del_len': {},
            'j_del_len': {},
            'vj_ins_len': {},
            'vdj_rarefaction': {}
        }
        __o_name = {}
        __sample_names = self.__sample_names
        __sample_names.append('all')
        for sample_name in __sample_names:
            sample_name_f = sample_name.replace(' ', '_')
            if sample_name not in __o_name:
                __o_name[sample_name] = {}
            for file_tag in __o_path.keys():
                _f = self.__output_prefix + '.' + sample_name_f + '.' + file_tag + '.freq.tsv'
                __o_path[file_tag][sample_name] = _f
                __o_name[sample_name][file_tag] = _f.split('/')[-1]
        self.__o_path = __o_path
        
        # Report
        self.__o_report = self.__output_prefix + '.report.html'
        
        # create objects
        stats = PyDAIRStats(self.__pydair_files, self.__sample_names,
                            self.__discard_ambiguous_D, self.__productive_only)
        
        if self.__args.estimate_vdj_combination:
            stats.rarefaction_study('vdj', self.__args.n_tries)
        
        self.stats  = stats
        self.report = PyDAIRReport(stats, __o_name)
    
    
    
    
    def write_summary(self, data_types = None):
        """Write summarized results.
        
        data_types (list): One of `v`, `d`, `j`, `vdj`, `cdr3_nucl`,
                           `cdr3_prot`, `v_del_len`, `j_del_len`,
                           and `vj_ins_len` can be specified.
                           If `None` is given, then write
                           all **data_types** results.
        
        Write the summarized into text file with TSV format.
        """
        
        if data_types is None:
            data_types = self.__o_path.keys()
        
        for data_type in data_types:
            for sample_i in range(len(self.stats.samples)):
                bsample = self.stats.samples.get_record(sample_i)
                data_type_freq = bsample.get_summary(data_type, func = 'raw')
                if data_type_freq is not None:
                    data_type_freq.to_csv(self.__o_path[data_type][self.__sample_names[sample_i]],
                                          sep = '\t', na_rep = 'NA',
                                          header = True, index = True)
                else:
                    with open(self.__o_path[data_type][self.__sample_names[sample_i]], 'w') as statsoutfh:
                        statsoutfh.write('Not avaliable')
            data_type_freq = self.stats.get_summary(data_type)
            
            if data_type_freq is not None:
                data_type_freq.to_csv(self.__o_path[data_type]['all'],
                                      sep = '\t', na_rep = 'NA',
                                      header = True, index = True)
            else:
                with open(self.__o_path[data_type][self.__sample_names[sample_i]], 'w') as statsoutfh:
                    statsoutfh.write('Not avaliable')
    
    
        
    
    def create_report(self):
        """Create report.
        
        Create an HTML report.
        """
        
        self.report.render(self.__o_report)
    
    
    
    

