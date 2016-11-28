======
Usages
======

Installation
============


PyDAIR can be installed via :command:`pip` command.

.. code-block:: bash

    pip install --user



Requirements:

* Python 2.7 or Python 3.4
* `NumPy <http://www.numpy.org/>`_
* `Pandas <http://pandas.pydata.org/>`_
* `BioPython <http://biopython.org/>`_
* `NCBI BLAST+ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_





parse action mode
=================

Usage
^^^^^

The required options in :command:`parse` action mode are 
``-q``, ``-v``, ``-d``, ``-j``, ``--v-blastdb``, ``--d-blastdb``, ``--j-blastdb``, and ``-o``.

* ``-q`` is used to specify input Ig-Seq data as FASTA format.
* ``-v``, ``-d``, and ``-j`` are used to specify FASTA format files of V, D, and J sequences.
* ``--v-blastdb``, ``--d-blastdb``, and ``--j-blastdb`` to specify BLAST database which created from the FASTA files by :command:`makeblastdb`.
* ``-o`` is used to specify the prefix of output files.


.. code-block:: bash
    
    pydair parse -q input_igh_sequences.fa \
                 -v v.fa                   \
                 -d d.fa                   \
                 -j j.fa                   \
                 --v-blastdb blastdb_v     \
                 --d-blastdb blastdb_d     \
                 --j-blastdb blastdb_j     \
                 -o output1



Options
^^^^^^^

In addition to the required options, 
the other options such as BLAST parameters for VDJ identification can be specified.


+----------------------------+------------------------------------+----------------+
| options                    |                                    | defult         |
+============================+====================================+================+
| ``-q``                     | Path to FASTA file that contains   |                |
|                            | IgH seqeunces.                     |                |
+----------------------------+------------------------------------+----------------+
| ``-o``                     | Path to file for writting results. |                |
+----------------------------+------------------------------------+----------------+
| ``-v``                     | Path to FASTA file of V gene       |                |
|                            | sequence.                          |                |
+----------------------------+------------------------------------+----------------+
| ``-d``                     | Path to FASTA file of D gene       |                |
|                            | sequence.                          |                |
+----------------------------+------------------------------------+----------------+
| ``-j``                     | Path to FASTA file of J gene       |                |
|                            | sequence.                          |                |
+----------------------------+------------------------------------+----------------+
| ``-v-blastdb``             | Path to BLAST database that        |                |
|                            | created from the FASTA file of     |                |
|                            | ``-v``.                            |                | 
+----------------------------+------------------------------------+----------------+
| ``-v-match-score``         | Score (> 0) for a nucleotide match | ``3``          |
|                            | for V gene.                        |                |
+----------------------------+------------------------------------+----------------+
| ``-v-mismatch-score``      | Score (< 0) for a nucleotide       | ``-3``         |
|                            | mismatch for V gene.               |                |
+----------------------------+------------------------------------+----------------+
| ``-v-gap-open-penalty``    | Penalty (> 0) to open a gap for    | ``6``          |
|                            | V gene.                            |                |
+----------------------------+------------------------------------+----------------+
| ``-v-gap-extend-penalty``  | Penalty (> 0) to extend a gap for  | ``6``          |
|                            | V gene.                            |                |
+----------------------------+------------------------------------+----------------+
| ``-v-wordsize``            | Word size to find hotspots by      | ``10``         |
|                            | BLAST for V gene.                  |                |
+----------------------------+------------------------------------+----------------+
| ``-v-evalue-cutoff``       | E-value                            | ``1e-10``      |
|                            | threshold for assiging V gene.     |                |   
+----------------------------+------------------------------------+----------------+
| ``-d-blastdb``             | Path to BLAST database that        |                |
|                            | created from the FASTA file of     |                |
|                            | ``-d``.                            |                | 
+----------------------------+------------------------------------+----------------+
| ``-d-match-score``         | Score (> 0) for a nucleotide match | ``1``          |
|                            | for D gene.                        |                |
+----------------------------+------------------------------------+----------------+
| ``-d-mismatch-score``      | Score (< 0) for a nucleotide       | ``-1``         |
|                            | mismatch for D gene.               |                |
+----------------------------+------------------------------------+----------------+
| ``-d-gap-open-penalty``    | Penalty (> 0) to open a gap for    | ``0``          |
|                            | D gene.                            |                |
+----------------------------+------------------------------------+----------------+
| ``-d-gap-extend-penalty``  | Penalty (> 0) to extend a gap for  | ``2``          |
|                            | D gene.                            |                |
+----------------------------+------------------------------------+----------------+
| ``-d-wordsize``            | Word size to find hotspots by      | ``4``          |
|                            | BLAST for D gene.                  |                |
+----------------------------+------------------------------------+----------------+
| ``-d-evalue-cutoff``       | E-value                            | ``1e-2``       |
|                            | threshold for assiging D gene.     |                |   
+----------------------------+------------------------------------+----------------+
| ``-j-blastdb``             | Path to BLAST database that        |                |
|                            | created from the FASTA file of     |                |
|                            | ``-j``.                            |                | 
+----------------------------+------------------------------------+----------------+
| ``-j-match-score``         | Score (> 0) for a nucleotide match | ``3``          |
|                            | for J gene.                        |                |
+----------------------------+------------------------------------+----------------+
| ``-j-mismatch-score``      | Score (< 0) for a nucleotide       | ``-3``         |
|                            | mismatch for J gene.               |                |
+----------------------------+------------------------------------+----------------+
| ``-j-gap-open-penalty``    | Penalty (> 0) to open a gap for    | ``6``          |
|                            | J gene.                            |                |
+----------------------------+------------------------------------+----------------+
| ``-j-gap-extend-penalty``  | Penalty (> 0) to extend a gap for  | ``6``          |
|                            | J gene.                            |                |
+----------------------------+------------------------------------+----------------+
| ``-j-wordsize``            | Word size to find hotspots by      | ``7``          |
|                            | BLAST for J gene.                  |                |
+----------------------------+------------------------------------+----------------+
| ``-j-evalue-cutoff``       | E-value                            | ``1e-5``       |
|                            | threshold for assiging J gene.     |                |   
+----------------------------+------------------------------------+----------------+
| ``-v-motif``               | The motif on V gene to determine   | ``YYC``        |
|                            | CDR3 segment.                      |                |   
+----------------------------+------------------------------------+----------------+
| ``-j-motif``               | The motif on J gene to determine   | ``WG.G``       |
|                            | CDR3 segment.                      |                |   
|                            | (e.g, ``WG.G``, ``FG.G``)          |                |   
+----------------------------+------------------------------------+----------------+



Output
^^^^^^


The :command:`parse` action generates some intermediate results with TSV format
and the final results with PYDAIR format.



+------------------------+------------+-----------------------------------------------------+
| file                   | format     | contents                                            |
+========================+============+=====================================================+
| <prefix>.v.blast.txt   | TSV        | BLAST results for V gene.                           |
+------------------------+------------+-----------------------------------------------------+
| <prefix>.j.blast.txt   | TSV        | BLAST results for J gene.                           |
+------------------------+------------+-----------------------------------------------------+
| <prefix>.vj.pydair     | PYDAIR     | Intermediate reuslts (V and J has been idenfied).   |
+------------------------+------------+-----------------------------------------------------+
| <prefix>.unaligned.fa  | FASTA      | Unaligned region sequences that will be used for    |
|                        |            | BLAST to identify D.                                |
+------------------------+------------+-----------------------------------------------------+
| <prefix>.d.blast.txt   | TSV        | BLAST results for D gene.                           |
+------------------------+------------+-----------------------------------------------------+
| <prefix>.vdj.pydair    | PYDAIR     | The final results of :command:`parse` mode.         |
+------------------------+------------+-----------------------------------------------------+




stats action mode
=================


Usage
^^^^^

Use ``-i`` to specify the PYDAIR format files generated by :command:`parse` action mode,
use ``-n`` to assigne the sample names to each PYDAIR format file,
and use ``-o`` to specify the prefix to save the summarized results.

.. code-block:: bash
    
    pydair stats -i output1.pydair output2.pydair output3.pydair  \
                 -n Fugu1 Fugu2 Fugu3                             \
                 -o stats_result                                  \
                 --contain_ambiguous_D



Option
^^^^^^

+--------------------------------+------------------------------------+----------------+
| options                        |                                    | defult         |
+================================+====================================+================+
| ``-i``                         | Path to ``PYDAIR`` format files.   |                |
|                                | Multiple files should be separated |                |
|                                | by a blank.                        |                |
+--------------------------------+------------------------------------+----------------+
| ``-n``                         | Sample names of each ``PYDAIR``    |                |
|                                | files. Multiple names should be    |                |
|                                | separated by a blank.              |                |
+--------------------------------+------------------------------------+----------------+
| ``-o``                         | Prefix for writting results.       |                |
+--------------------------------+------------------------------------+----------------+
| ``--contain_ambiguous_D``      | If ``True``, summarize all         | ``True``       |
|                                | sequences regardless the D segment |                |
|                                | is identified or not.              |                |
|                                | If ``False``, summarize only the   |                |
|                                | seqeunces with identified D        |                |
|                                | segment.                           |                |
+--------------------------------+------------------------------------+----------------+
| ``--productive_only``          | If ``False``, summarize all        | ``False``      |
|                                | sequences regardless productive or |                |
|                                | nonproductive sequences.           |                |
|                                | If ``True``, summarize only the    |                |
|                                | productive sequences.              |                |
+--------------------------------+------------------------------------+----------------+
| ``--estimate-vdj-combination`` | If ``True``, perform rarefaction   | ``False``      |
|                                | analysis to study the diversity of |                |
|                                | VDJ combinations.                  |                |
+--------------------------------+------------------------------------+----------------+


Note that, the productive sequence is defined as the sequence without any stop codons,
whereas the nonproductive sequence is defined as the sequence with at least one stop codons.



Output
^^^^^^


+-------------------------------------+------------+-----------------------------------------------------+
| file                                | format     | contents                                            |
+=====================================+============+=====================================================+
| <prefix>.v.freq.tsv                 | TSV        | V gene usage frequency.                             |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.d.freq.tsv                 | TSV        | D gene usage frequency.                             |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.j.freq.tsv                 | TSV        | J gene usage frequency.                             |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.vdj.freq.tsv               | TSV        | Frequencies of VDJ combinations.                    |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.cdr3_nucl_length.freq.tsv  | TSV        | CDR3 nucleotide sequence distribution.              |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.cdr3_prot_length.freq.tsv  | TSV        | CDR3 amino acid sequence distribution.              |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.rarefaction.tsv            | TSV        | Rarefaction analysis results of VDJ combination.    |
+-------------------------------------+------------+-----------------------------------------------------+
| <prefix>.report.html                | HMLT       | HTML report of summarized results.                  |
+-------------------------------------+------------+-----------------------------------------------------+






PYDAIR format
=============

PYDAIR format is an human readable text file format.
Each entry in PYDAIR format file represents a result of
VDJ identification and CDR3 segment identification for a query IgH sequence. 
PYDAIR format file can contain multiple entries.
In PYDAIR format, an entry is begin with ``#BEGIN`` and finished with ``#END``.


.. code-block:: text
    
    #BEGIN
    QN M03400:8:000000000-ADYYJ:1:1101:12827:2301
    VN v1.14
    DN .
    JN Jm1
    OP 1
    OC _M____;**__**;**__*M
    QA                  CTGACCCAGTCTGAACCAGTGGTTAAAAGACCTGGAGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGCAGCTATGGGATGAACTGGGTCAGACAGGCTCCTGGAAAAGGACTGGAGTGGATTGCTTATATCTAC------AGCAGCACATACTACTCTGAGTCAGTCAAAGGCCGGTTTAACATCTCCAGAGATAACAACAGAGCACAGCTGAATCTGCATATAAAAAGCCTGAAGACTGAAGATACTGCGGGTTATTATTGTGCTCGAACTGGAAAAGAATACCTTGACTACTGGGGGAAAGGCACAACAGTTACAGTAACGTCTGCAACACCCAAAGCCCCTTCTTGGTTTCCATTGATACAATGCGGAACTGGGACTGGAACCCTGGTCACTCTCGGATGTTTGGCCGCCGACTTCACGCCATCGGACCTAACCTACACCTGGAGAAAAGACGGAGTCGATCTGAAAGACTTCATTCAGTACCCTCCAACCACGAACG
    VA GTGTTGATGCTCAGACTCTGACCCAGTCTGAACCAGTGGTTAAAAGACCTGGGGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGCAGCTACTGGATGGTCTGGGTCAGACAGGCTCCTGGAAAAGGACTGGAGTGGATCGCTTATATCACCACCAGTAGCAGCCCATACTACTCTGAGTCAGTCAAAGGCCGGTTTATCATCTCCAGAGACAACAACAGAGCACAGCTGAATCTGCAGATTAACAGCCTGAAGACTGAAGATTCTGCTGTTTATTATTGTGCTCGAGAG
    JA                                                                                                                                                                                                                                                                                                               TACTACGCATACTTTGACTACTGGGGGAAAGGAACAACAGTTACAGTAACATCT
    UA                                                                                                                                                                                                                                                                                                             CTGGAAAAGA
    CA                                                                                                                                                                                                                                                                                                      GCTCGAACTGGAAAAGAATACCTTGACTACTGG
    #CDR3AA ARTGKEYLDYW
    #AL QSTART  QEND    SSTART  SEND    IDENTITY    SCORE
    AL QV 1 276 18  299 91.844  1.29e-100
    AL QD . .   .   .   .   .
    AL QJ 288   333 9   54  93.478  5.62e-15
    AL QU 277   287 .   .   .   .
    AL QC 270   303 .   .   .   .
    #END



Each line begins with a 2-character line code
that indicates the type of information contained in the line.
There are 13 line codes defined.
In addition, the line begining with a ``#`` is a comment line,
which gives additional information but not required.


+-----------+------------------------------------------------+
| line code | definition                                     |
+===========+================================================+
| ``QN``    | Query sequence name.                           |
+-----------+------------------------------------------------+
| ``VN``    | Assigned V gene name.                          |
+-----------+------------------------------------------------+
| ``DN``    | Assigned D gene name.                          |
+-----------+------------------------------------------------+
| ``JN``    | Assigned J gene name.                          |
+-----------+------------------------------------------------+
| ``OP``    | The start position of reading frame of         |
|           | the query sequence.                            |
+-----------+------------------------------------------------+
| ``OC``    | The information about start codons and         |
|           | stop codons.                                   |
+-----------+------------------------------------------------+
| ``QA``    | The aligned sequence of query.                 |
+-----------+------------------------------------------------+
| ``VA``    | The aligned sequence of V gene.                |
+-----------+------------------------------------------------+
| ``JA``    | The aligned sequence of J gene.                |
+-----------+------------------------------------------------+
| ``UA``    | The aligned sequence of un-aligned region.     |
+-----------+------------------------------------------------+
| ``CA``    | The aligned sequence of CDR3.                  |
+-----------+------------------------------------------------+
| ``AL``    | The inforamtion about BLAST results and        |
|           | alignment information. This line code          |
|           | consists of five attributes.                   |
+-----------+------------------------------------------------+

There are five attributes defined in ``AL``.

+------------+--------------------------------------------------+
| attribute  | definition                                       |
+============+==================================================+
| ``QV``     | Alignment positions and BLAST statistics between |
|            | a query and the V gene.                          |
+------------+--------------------------------------------------+
| ``QD``     | Alignment positions and BLAST statistics between |
|            | a query and the D gene.                          |
+------------+--------------------------------------------------+
| ``QJ``     | Alignment positions and BLAST statistics between |
|            | a query and the J gene.                          |
+------------+--------------------------------------------------+
| ``QU``     | Alignment positions between a query and the      |
|            | unaligned region.                                |
+------------+--------------------------------------------------+
| ``QC``     | Alignment positions between a query and the      |
|            | CDR3 sequence.                                   |
+------------+--------------------------------------------------+

Each attribute contains the five columns that separated by TAB.
The five columns indicates that the alignement start and end
positions in query sequence,
the alignemnt start and end positions in subject (V, D, J, unaligned, and CDR3) sequence,
and the identity and BLAST score of alignemnt.







