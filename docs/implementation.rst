==============
Implementation
==============



parse mode
==========

The parse mode is designed for identifying VDJ gene segments and CDR-H3 region,
and detecting indel among VDJ junctions. PyDAIR internally calls BLAST.
User-defined VDJ database are required prior to the use of BLAST. 

.. image:: images/pydair-parse.png
    :align: center







#. **V and J assignment** 
   All sequences were aligned by using BLAST (i.e, blastn on NCBI BLAST+)
   to use-defined V- and J- gene segment databases with match gain of +3,
   mismatch cost of -3, gap-open penalty of 6, and gap-extend penalty of 6.
   The optimal alignments for V- (V-aligned region) and J- (J-aligned region)
   gene segments were assigned for each IgH sequence.
   Only sequences with determined V- and –J assignments were kept for downstream analysis. 
 
#. **CDR-H3 determination**
   CDR-H3 is defined as a region between conserved residues of V (C in 2nd-CYS motif)
   and J (W in WGxG motif of J) (North et al. 2011; Giudicelli 2006; Savan et al. 2005).
   Each IgH sequence was translated into amino acid (aa) sequence,
   and was searched for the two motifs using regular expression matching.
   In detail, the YYC motif was sought from the 5’ end of the J-aligned region
   in 3’ to 5’ direction and the WGxG motif was searched from the 3’ end of the
   V-aligned region in 5’ to 3’ direction, as illustrated in Figure S3 B.
   It is worth noting that the motif searches were executed separately,
   that iteration was performed for the three possible reading frames one at a time.
   It is therefore possible that the reading frame of YYC differs from that of WGxG,
   which results in out-of-frame CDR-H3 nucleotide (nt) sequence (length of sequence
   is not a multiple of three).
   In this case, CDR-H3 aa sequence was defined as the sequence translating from the
   first nucleotide from CDR-H3 nt sequence.
   In addition, PyDAIR allows users to define species-specific motifs when determining
   CDR3 region by adjusting ‘--v-motif’ and ‘--j-motif’ options.
    .. image:: images/pydair-parse-motif.png
        :align: center

#. **Indel detection**
   As shown in Figure S3 C, once the V- or J- aligned region has been assigned
   (the local alignment), PyDAIR adds expected nucleotides of V or J gene sequence
   adjacent to the local alignment to form a global alignment.
   The number of deletions is calculated by counting bases starting from the end of
   the last three consecutive nucleotides within global alignment and finishing at
   the last nucleotide of expected V or J gene.
   Finally, the sequence found between the 3’ end of V and the 5’ end of J is defined
   as a non-template V-J insertion region. 
    
    .. image:: images/pydair-parse-indels.png
        :align: center

#. **D assignment**
   The V-J insertion region contains the D gene segment used.
   PyDAIR employs the V-J insertion region for aligning with the user-defined D gene database.
   A V-J insertion sequence was retained if it contained at least 4 nucleotides
   for D gene identification. Parameters used were shown in Table S2. 

#. **Reading frame detection**
   Open reading frame (ORF) sequence was determined for each alignment by iteratively
   searching for stop codons (i.e., TAG, TAA, and TGA) on the three possible reading
   frames with regular expression matching.
   Iteration was performed until at least one of the reading frames containing no stop
   codons could be found.
   At this point, the IgH sequence was defined as productive, otherwise,
   ‘unproductive’ was assigned. 



stats mode
==========

The stats mode provides a number of functions for summarizing and visualizing analysis
results from pars mode; it summarizes:
(i) usage frequencies of VDJ gene segments;
(ii) frequencies of VDJ combinations, rarefaction analysis for estimating the
saturation of VDJ combinations; 
(iii) length distributions of CDR-H3 nt and aa sequences;
and (iv) indel distribution within V-D and D-J junctions.
Finally, the mode generates an HTML report by incorporating
plot.ly API (https://plot.ly/) for visualizing summarized results.



sim mode
========

The sim mode is used for artificial IgH repertoire generation.
In some cases, BLAST might give weak results when the local database consisting
of highly similar sequences, thus,
it is recommended to perform procedures in sim and eval modes for evaluating
performance of PyDAIR prior to data analysis. 

The sim mode generates artificial IgH repertoire as follows:

1) Randomly sample each V, D, and J gene segment from user-defined databases.
2) Trim bases at each putative junction boundary i.e., 5’- and 3’- end of
   sampled V, D, and J gene segments. Default base number for trimming is shown in Table S5. 
3) Randomly insert sampled bases from A, T, C, G into VDJ junctions wherein
   VD represents bases being added to V-D junction and aDJ to D-J junction (Table S5). 
4) Mutate each nucleotide in the sequence with probability of p (default 0.05).
   Probabilities of mutation (substitution or deletion) are set to equal frequency
   (i.e., the probability for each mutation is of 0.25).


eval mode
=========

Evaluation of PyDAIR performance is performed using the eval mode by calculating
the number of (in)correctly IgH sequences assigned for certain VDJ gene segments. 




    



