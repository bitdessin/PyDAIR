==============
Implementation
==============


There four main modes (i.e., functions) are implemented in PyDAIR:
:command:`parse`, :command:`stats`, :command:`sim`, and :command:`eval` modes.


+---------------------+------------------------------------------------------+
| Mode                | Functions                                            |
+=====================+======================================================+
| :command:`parse`    | V, D, J, and CDR3 segments identification and        |
|                     | indels detection among VD and DJ junctions.          |
+---------------------+------------------------------------------------------+
| :command:`stats`    | Results summarizaiton and visualization.             |
|                     | This mode should be used after :command:`parse` mode.|
+---------------------+------------------------------------------------------+
| :command:`sim`      | Generating artificial IgH sequences for simulaiton   |
|                     | studies.                                             |
+---------------------+------------------------------------------------------+
| :command:`eval`     | Evaluating PyDAIR performances with artificial       |
|                     | IgH sequencs. This mode should be used after         |
|                     | :command:`sim` and :command:`parse` modes.           |
+---------------------+------------------------------------------------------+




parse mode
==========

The :command:`parse` mode is designed for identifying V, D, J, and CDR3 segments,
and detecting indels among VD and DJ junctions.
NCBI BLAST+ (:command:`blastn`) is internally used in PyDAIR :command:`parse` mode.


.. image:: images/pydair-parse.png
    :align: center


#. **V and J assignment** 
   All sequences were aligned by using BLAST (:command:`blastn`)
   to user-defined V and J gene segment databases.
   The optimal alignments for V (V-aligned region) and J (J-aligned region)
   gene segments were assigned for each IgH sequence.
   Only sequences with determined V and J assignments were kept for downstream analysis.
   PyDAIR allows users to change BLAST parameters (see :doc:`Usages<usage>`).
   
 
#. **CDR3 determination**
   CDR3 is defined as a region between the conserved residues of V (YYC motif) and
   J (WGxG motif). CDR3 is determined as the sequences between the motifs.
   In detail, the YYC motif is sought from the 5' end of the J-aligned region
   in the 3' to 5' direction,
   and the WGxG motif is searched from the 3' end of the V-aligned region
   in the 5' to 3' direction.
   Notably, the motif searches were executed separately,
   and iteration was performed for three possible reading frames, one at a time.
   It is therefore possible that the reading frame of YYC differs from that of WGxG,
   resulting in out-of-frame CDR3 nucleotide (nt) sequences
   (i.e., the sequence length is not a multiple of three).
   In this case, the CDR3 aa sequence was defined as the sequence translated from
   the first nucleotide of the CDR3 nt sequence.
   PyDAIR allows users to define species-specific motifs when determining
   CDR3 region (see :doc:`Usages<usage>`).
    
    .. image:: images/pydair-parse-motif.png
        :align: center

   
#. **Indel detection** 
   Once the V- or J-aligned region has been assigned (*local alignment*),
   PyDAIR adds expected nucleotides of V or J gene sequence
   adjacent to the local alignment to form a *global alignment*.
   The number of deletions is calculated by counting bases starting from the end of
   the last three consecutive nucleotides within global alignment and finishing at
   the last nucleotide of expected V or J gene.
   Finally, the sequence found between the 3' end of V and the 5' end of J is defined
   as a non-template V-J insertion region. 
    
    .. image:: images/pydair-parse-indels.png
        :align: center

   
#. **D assignment**
   The V-J insertion region, with at least 4 nucelotides, is used for
   assignment of D gene segment.


#. **Reading frame detection**
   Stop codons (i.e., TAG, TAA, and TGA) are searched from whole IgH sequence with
   the three reading frames.
   If at least one frame dose not contain stop codon, the sequence
   is defiend as *productive*, otherwise, defined as *uproductive*.
    




stats mode
==========

The :command:`stats` mode provides a number of functions for summarizing and visualizing 
the analysis results from :command:`parse` mode; it summarizes:
(i) usage frequencies of VDJ gene segments;
(ii) frequencies of VDJ combinations including rarefaction analysis for estimating the
saturation of VDJ combinations; 
(iii) length distributions of CDR3 nt and aa sequences;
and (iv) indel distributions among V-D and D-J junctions.
The summarization are generated as HTML report using `plot.ly <https://plot.ly/>`_.



sim mode
========

Since BLAST might have inaccuracy results
when the local databases consisting of highly similar sequences,
it is recommended to evaluate the performance of PyDAIR with artificial sequences
prior to data analysis.
There are two commands,  :command:`sim` and :command:`eval`, used for this purpose.

The :command:`sim` mode is used for generating artificial IgH sequences
with the following steps.

1) Randomly sample each V, D, and J gene segment from user-specified sequences.
2) Trim several nucleotides at each putative junction boundary
   i.e., 5' and 3' end of sampled V, D, and J gene segments.
   The number of nucleotides for trimming can be specified with parameters.
3) Randomly insert sampled nucleotides from A, T, C, G into V-D and D-J junctions.
   The number of nucleotides for insertion can be specified with parameters.
4) Mutate each nucleotide in the sequence with probability of *p*.
   Probabilities of mutation (substitution or deletion) are set to equal frequency
   (i.e., the probability for each mutation is of 0.25).


eval mode
=========

Evaluation of PyDAIR performances is performed using the :command:`eval` mode by calculating
the number of correctly and incorrectly IgH sequence assignments.




