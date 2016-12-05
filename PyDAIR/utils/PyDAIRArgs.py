import os
import argparse





class PyDAIRParseSeqArgs:
    """A class to save parameters for parsing Rep-Seq data.
    
    PyDAIRParseSeqArgs is a class to store the parameters for parsing Rep-Seq data.
    """
   
    def __init__(self, q = None, v = None, d = None, j = None, o = None,
                 valign = None, dalign = None, jalign = None,
                 v_motif = None, j_motif = None):
        """Set up parameters for analysis of Rep-Seq data
        
        Args:
            q (str): A path to V gene FASTA file.
            d (str): A path to D gene FASTA file.
            j (str): A path to J gene FASTA file.
            o (str): A prefix to save the result files.
            valign (PyDAIRBlastArgs): BLAST parameters for V gene identification.
            dalign (PyDAIRBlastArgs): BLAST parameters for D gene identification.
            jalign (PyDAIRBlastArgs): BLAST parameters for J gene identification.
            v_motif (str): The motif on V gene to determine CDR3 segment.
            j_motif (str): The motif on J gene to determine CDR3 segment.
        
        Raises:
            IOError: An error occurs when the FASTA files are not existed.
        
        >>> valign = PyDAIRBlastArgs('balstdb_v', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 21, eval_cutoff = 1e10)
        >>> dalign = PyDAIRBlastArgs('balstdb_d', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 2, eval_cutoff = 1e3)
        >>> jalign = PyDAIRBlastArgs('balstdb_j', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 9, eval_cutoff = 1e10)
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
    """A class to save BLAST parameters for V, D, J gene identification.
    
    PyDAIRBlastArgs is a class to save BLAST parameters for V, D, J gene 
    identification.
    Usually, the three PyDAIRBlastArgs class objects (for V, D, and J genes)
    are used simulatory.
    """
    
    def __init__(self, db, match, mismatch, gapopen, gapextend, wordsize, eval_cutoff):
        """Set up parameters for statistical analysis.
        
        Args:
            db (str): A path to BLAST database name. 
                      The BLAST database is made by ``makeblastdb`` command in BLAST program.
            match (int): A positive integer represents the score for a nucleotide match.
                         The argument will pass to ``-reward`` in BLAST.
            mismatch (int): A negative integer represents the score for a nucleotide match.
                            The argument will pass to ``-penalty`` in BLAST.
            gapopen (int): A positive integer represents the penalty to open a gap.
                           The argument will pass to ``-gapopen`` in BLAST.
            gapextend (int): A positive integer represents the penalty to extend a gap.
                             The argument will pass to ``-gapopen`` in BLAST.
            wordsize (int): A positive integer represents word size.
                            The argument will pass to ``-word_size`` in BLAST.
            eval_cutoff (float): The cutoff of e-value for identifying genes.
        
        >>> valign = PyDAIRBlastArgs('balstdb_v', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 21, eval_cutoff = 1e10)
        >>> dalign = PyDAIRBlastArgs('balstdb_d', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 2, eval_cutoff = 1e3)
        >>> jalign = PyDAIRBlastArgs('balstdb_j', match = 5, mismatch = 6, gapopen = 4,
        >>>                          gapextend = 4, wordsize = 9, eval_cutoff = 1e10)
        
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
    """A class is to save parameters for statistical analysis.
    
    There two arguments `discard_ambiguous_D` and `productive_only`
    are important in analysis.
    
    An Ig sequence is composed of V, D, and J genes. In addition,
    The region between the end of V (known as YYC motif) and
    the start of J (known as WGxG motif) is defined CDR3.
    During the VDJ recombination processes,
    number of additions, deletions, and mutations occurs on CDR3.
    Therefore, it is difficult to identify D gene usage in Ig seqeunce.
    
    After parsing (identification of V, D, J genes and determiniation of
    CDR3 seqeunce) Rep-Seq data, there plenty of Ig sequence may have
    unidentificable D genes (i.e., ambiguous D).
    
    The `discard_ambiguous_D` argument is an option to specify
    how to treat these Ig sequences with ambiguous D gene.
    If `discard_ambiguous_D` is `True`, analysis will be performed
    against all Ig sequences.
    If `discard_ambiguous_D` is `False`, before analysis,
    all Ig sequences with ambiguous D gene are discarded.
    Then, analysis will be performed agains the remained Ig
    sequences with identificable D gene.

    The `productive_only` argument is an option to specify
    how to treat these Ig seqeunces that contains stop codons.
    It is known to that Ig seqeunce freaquently contains stop codons.
    The Ig sequence with stop codons may not give the function in immune system.
    If `productive_only` is `True`, analysis will be performed against
    all Ig sequences.
    If `productive_only` is `False`, analysis will be performed against
    Ig seqeuences that do not contain any stop codon.
    """
    
    def __init__(self, sample_names, pydair_files, discard_ambiguous_D, productive_only,
                 estimate_vdj_combination, n_tries = 1000,
                 output_prefix = './pydairstats_'):
        """Set up parameters for statistical analysis.

        Args:
            sample_names (list): A list of string represents the sample names.
                                 The analysis results will use ``sample_names`` instead of file names.
            pydair_files (list): A list of string containing path to PYDAIR flat file.
                                 The length of ``pydair_files`` should be equal to ``sample_names``.
            discard_ambiguous_D (bool):  If ``True``, analysis will contain Ig sequence with ambiguous D gene.
                                         Default is ``True``.
            productive_only (bool): If ``True``, analysis will contain Ig sequence which contains stop codons.
                                      Since Ig sequence with stop codon may not have function in immune system,
                                      therefore, default is ``False``.
            estimate_vdj_combination (bool): If ``True``, perform rarefaction study to estimate diversity of VDJ combination.
            output_prefix (str): A prefix to save the analysis results.
            n_tries (int): Number of simulation tries.
        
        >>> args = PyDAIRStatsArgs(sample_names = ['fugu 1', 'fugu 2', 'fugu 3'],
        >>>                        pydair_files = ['fugu1.pydair', 'fugu2.pydair', 'fugu3.pydair'],
        >>>                        containe_ambiguous_D = True, productive_only = False,
        >>>                        estimate_vdj_combination = True, n_tries = 1000,
        >>>                        output_prefix = './fugustat_result')
        
        """
        
        self.sample_names        = sample_names
        self.pydair_files        = pydair_files
        self.discard_ambiguous_D = discard_ambiguous_D
        self.productive_only   = productive_only
        self.output_prefix       = output_prefix
        #self.figure_format       = figure_format
        #self.figure_style        = figure_style
        self.estimate_vdj_combination = estimate_vdj_combination
        self.n_tries             = n_tries

