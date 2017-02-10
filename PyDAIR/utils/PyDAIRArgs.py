import os
import argparse





class PyDAIRParseSeqArgs:
    """PyDAIRParseSeqArgs class.

    PyDAIRParseSeqArgs class is used for saving the parameters for
    parsing Rep-Seq data.
    """
   
    def __init__(self, q = None, v = None, d = None, j = None, o = None,
                 valign = None, dalign = None, jalign = None,
                 v_motif = None, j_motif = None):
        """Set up parameters for parsing Rep-Seq data.
        
        Args:
            q (str): A path to V gene FASTA file.
            d (str): A path to D gene FASTA file.
            j (str): A path to J gene FASTA file.
            o (str): A prefix to save the result files.
            valign (PyDAIRBlastArgs): BLAST parameters for V segment identification.
            dalign (PyDAIRBlastArgs): BLAST parameters for D segment identification.
            jalign (PyDAIRBlastArgs): BLAST parameters for J segment identification.
            v_motif (str): The motif on 3'-end of V gene for CDR3 segment identification.
            j_motif (str): The motif on 5'-end of J gene for CDR3 segment identification.
        
        
        >>> valign = PyDAIRBlastArgs('balstdb_v', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 21, eval_cutoff = 1e-10)
        >>> dalign = PyDAIRBlastArgs('balstdb_d', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 2, eval_cutoff = 1e-3)
        >>> jalign = PyDAIRBlastArgs('balstdb_j', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 9, eval_cutoff = 1e-10)
        >>> args = PyDAIRParseSeqArgs('fugu', q = 'q.fa', v = 'v.fa', d = 'd.fa',
        >>>                           j = 'j.fa', o = 'out_prefix', f = 'pydair',
        >>>                           valign = valign, dalign = dalign, jalign = jalign,
        >>>                           v_motif = 'YYC', j_motif = 'WG.G')
        
        """

        if os.path.exists(q) is False:
            raise IOError('File [' + q + '] not found.')
        if os.path.exists(v) is False:
            raise IOError('File [' + v + '] not found.')
        if os.path.exists(d) is False:
            raise IOError('File [' + d + '] not found.')
        if os.path.exists(j) is False:
            raise IOError('File [' + j + '] not found.')
        
        self.q_file_path = q
        self.v_file_path = v
        self.d_file_path = d
        self.j_file_path = j
        self.output_path = o
        self.v_align_args = valign
        self.d_align_args = dalign
        self.j_align_args = jalign
        self.v_motif = v_motif
        self.j_motif = j_motif



class PyDAIRBlastArgs:
    """PyDAIRBlastArgs class.
    
    PyDAIRBlastArgs class is used for saving the BLAST parameters for
    V, D, and J segment identification.
    
    """
    
    def __init__(self, db, match, mismatch, gapopen, gapextend, wordsize, eval_cutoff):
        """Set up BLAST parameters for V, D, and J identification.
        
        Args:
            db (str): A path to BLAST database which created by `makeblastdb`. 
            match (int): A positive integer represents the score for a nucleotide match.
                         The value will pass to `-reward` in BLAST.
            mismatch (int): A negative integer represents the score for a nucleotide match.
                            The value will pass to `-penalty` in BLAST.
            gapopen (int): A positive integer represents the penalty to open a gap.
                           The value will pass to `-gapopen` in BLAST.
            gapextend (int): A positive integer represents the penalty to extend a gap.
                             The value will pass to `-gapopen` in BLAST.
            wordsize (int): A positive integer represents word size.
                            The value will pass to `-word_size` in BLAST.
            eval_cutoff (float): A decimal used as E-value cutoff for
                                 segment identificaiton.

        >>> valign = PyDAIRBlastArgs('balstdb_v', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 21, eval_cutoff = 1e-10)
        >>> dalign = PyDAIRBlastArgs('balstdb_d', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 2, eval_cutoff = 1e-3)
        >>> jalign = PyDAIRBlastArgs('balstdb_j', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 9, eval_cutoff = 1e-10)
        
        """
        
        if match < 0:
            raise argparse.ArgumentTypeError('The match-score argument for V, D, and J should be a positive integer.')
        if mismatch > 0:
            raise argparse.ArgumentTypeError('The mismatch-score argument for V, D, and J should be a negative integer.')
        if gapopen < 0:
            raise argparse.ArgumentTypeError('The gap-open-penalty argument for V, D, and J should be a positive integer.')
        if gapextend < 0:
            raise argparse.ArgumentTypeError('The gap-extend-penalty argument for V, D, and J should be a positive integer.')
        if wordsize < 3:
            raise argparse.ArgumentTypeError('The wordsize argument for V, D, and J should be positive integer that is not less than 4.')
        self.db        = db
        self.match     = str(match)
        self.mismatch  = str(mismatch)
        self.gapopen   = str(gapopen)
        self.gapextend = str(gapextend)
        self.wordsize  = str(wordsize)
        self.cutoff    = float(eval_cutoff)





class PyDAIRStatsArgs:
    """PyDAIRStatsArgs class.
    
    PyDAIRStatsArgs class is used for saving parameters for summarizing parsed results.
    
    """
    
    def __init__(self, sample_names, pydair_files,
                 discard_ambiguous_D = False,
                 productive_only = False,
                 estimate_vdj_combination = False,
                 n_tries = 1000,
                 output_prefix = './pydairstats_'):
        """Set up parameters for summarizing parsed results.

        Args:
            sample_names (list): A list of string represents the sample names.
            pydair_files (list): A list of string containing path to PYDAIR flat file.
                                 The length of **pydair_files** should equal to
                                 that of **sample_names**.
            discard_ambiguous_D (bool): Default is `False`. If `True`, the sequences
                                        with ambiguous D segment will be discarded before
                                        summarization.
            productive_only (bool): Default is `False`. If `True`, the sequences contained
                                    more than one stop codon will be discarded before
                                    summarization.
            estimate_vdj_combination (bool): Default is `False`. If `True`, performe
                                             rarefaction analysis to estimate the diversity
                                             of VDJ combination.
            n_tries (int): Number of simulation tries for rarefaction analysis.
            output_prefix (str): A prefix to save the result files.
        
        >>> args = PyDAIRStatsArgs(sample_names = ['fugu 1', 'fugu 2', 'fugu 3'],
        >>>                        pydair_files = ['fugu1.pydair', 'fugu2.pydair', 'fugu3.pydair'],
        >>>                        containe_ambiguous_D = True,
        >>>                        productive_only = True,
        >>>                        estimate_vdj_combination = True, n_tries = 1000,
        >>>                        output_prefix = './fugustat_result')
        
        """
        
        self.sample_names        = sample_names
        self.pydair_files        = pydair_files
        self.discard_ambiguous_D = discard_ambiguous_D
        self.productive_only   = productive_only
        self.output_prefix       = output_prefix
        self.estimate_vdj_combination = estimate_vdj_combination
        self.n_tries             = n_tries




class PyDAIRSimArgs:
    """PyDAIRSimArgs class.
    
    PyDAIRSimArgs class is used for saving parameters for generation of artificial
    IgH sequences.
    
    """
    
    def __init__(self, o, n,
                       v_fa, n_v_5del, n_v_3del,
                       d_fa, n_d_5del, n_d_3del,
                       j_fa, n_j_5del, n_j_3del,
                       n_vd_ins, n_dj_ins, p_mutation, random_seed = None):
        """Set up parameters for generation of artificial IgH sequences.
        
        """
        
        self.output = o
        self.n = n
        self.v_fasta = v_fa
        self.d_fasta = d_fa
        self.j_fasta = j_fa
        self.n_v_5del = n_v_5del
        self.n_v_3del = n_v_3del
        self.n_d_5del = n_d_5del
        self.n_d_3del = n_d_3del
        self.n_j_5del = n_j_5del
        self.n_j_3del = n_j_3del
        self.n_vd_ins = n_vd_ins
        self.n_dj_ins = n_dj_ins
        self.p_mutation = p_mutation
        self.seed = random_seed




class PyDAIREvalArgs:
    """PyDAIREvalArgs class.
    
    """
    
    def __init__(self, sim_condition, parse_result, output):
        self.output = output
        self.sim_condition = sim_condition
        self.parse_result = parse_result

