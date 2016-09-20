import os
import argparse





class PyDAIRParseSeqArgs:
    '''
    PyDAIRParseSeqArgs class is a class to store the parameters for analyzing data.
    The initialize method will check that the file path of 'q', 'v', 'd', and 'j'. 
    
      q:        FASTA file of query sequences.
      v:        FASTA file of V gene sequences.
      d:        FASTA file of D gene sequences.
      j:        FASTA file of J gene sequences.
      o:        The prefix of output file.
      f:        PyDAIR format, 'PyDAIR', 'tsv', and 'simple' can be specified.
      valign:   PyDAIRBlastArgs class object for V gene.
      dalign:   PyDAIRBlastArgs class object for V gene.
      jalign:   PyDAIRBlastArgs class object for V gene.
    
    '''
    def __init__(self, species = None, q = None, v = None, d = None, j = None, o = None, f = None,
                 valign = None, dalign = None, jalign = None):
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
        self.pydair_format = f
        self.v_align_args = valign
        self.d_align_args = dalign
        self.j_align_args = jalign
        self.species = species




class PyDAIRBlastArgs:
    '''
    PyDAIRBlastArgs is a class to store the parameters for BLAST research. Usually, the three
    PyDAIRBlastArgs (for V, D, and J genes) should be used simulatory.
    The initialization method will check the several parameters
    as metioned as follwoings.
    
      db:          The path to the BLAST database which made by 'makeblastdb' command in BLAST program.
      match:       Match score. Positive integer.
      mismatch:    Mismatch score. Negative integer.
      gapopen:     Gap-open penalty. Negative integer.
      gapextend:   Gap-extend penalty. Negative integer.
      wordsize:    Wordsize for BLAST. Positive integer, and it should be larger than 3.
      eval_cutoff: The cutoff of e-value for identifying genes.
    
    '''
    def __init__(self, db, match, mismatch, gapopen, gapextend, wordsize, eval_cutoff):
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
    '''
    PyDAIRStatsArgs is a class to store parameters for statistical analyzing.
    
      sample_names: Sample 
    
    '''
    def __init__(self, sample_names, pydair_files, contain_ambiguous_D, contain_stopcodon,
                 output_prefix, figure_format = 'pdf', figure_style = 'ggplot'):
        self.sample_names        = sample_names
        self.pydair_files         = pydair_files
        self.contain_ambiguous_D = contain_ambiguous_D
        self.contain_stopcodon   = contain_stopcodon
        self.output_prefix       = output_prefix
        self.figure_format       = figure_format
        self.figure_style        = figure_style
        










