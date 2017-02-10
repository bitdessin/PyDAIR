======
Usages
======

Installation
============


PyDAIR is a Python package which can be installed via :command:`pip` command.

.. code-block:: bash

    pip install pydair --user


Requirements:

* Python 2.7 or Python 3.4
* `NumPy <http://www.numpy.org/>`_ (Python package)
* `Pandas <http://pandas.pydata.org/>`_ (Python package)
* `BioPython <http://biopython.org/>`_ (Python package)
* `NCBI BLAST+ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_





parse action mode
=================

The :command:`parse` action mode identifies V, D, J, and CDR3 segments,
and detects the deletions and insertions in IgH sequences.

Usage
^^^^^

.. code-block:: bash
    
    pydair parse -q input_igh_sequences.fa \
                 -o output1                \
                 -v v.fa                   \
                 -d d.fa                   \
                 -j j.fa                   \
                 --v-blastdb blastdb_v     \
                 --d-blastdb blastdb_d     \
                 --j-blastdb blastdb_j 


The eight options, i.e., ``-q``, ``-o``, ``-v``, ``-d``, ``-j``,
``--v-blastdb``, ``--d-blastdb``, and ``--j-blastdb``, are required options.


Options
^^^^^^^

In addition to the required options, 
the other options such as BLAST parameters for VDJ identification can be specified.
All options used in :command:`parse` action mode are shown in the following table.


+-----------------------------+------------------------------------+----------------+
| Options                     | Description                        | Default        |
+=============================+====================================+================+
| ``-q``                      | Path to FASTA file that contains   |                |
|                             | IgH seqeunces.                     |                |
+-----------------------------+------------------------------------+----------------+
| ``-o``                      | Path to file for writting results. |                |
+-----------------------------+------------------------------------+----------------+
| ``-v``                      | Path to FASTA file of V gene       |                |
|                             | sequence.                          |                |
+-----------------------------+------------------------------------+----------------+
| ``-d``                      | Path to FASTA file of D gene       |                |
|                             | sequence.                          |                |
+-----------------------------+------------------------------------+----------------+
| ``-j``                      | Path to FASTA file of J gene       |                |
|                             | sequence.                          |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-blastdb``             | Path to BLAST database that        |                |
|                             | created from the FASTA file of     |                |
|                             | ``-v``.                            |                | 
+-----------------------------+------------------------------------+----------------+
| ``--v-match-score``         | Score (> 0) for a nucleotide match | ``3``          |
|                             | for V gene.                        |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-mismatch-score``      | Score (< 0) for a nucleotide       | ``-3``         |
|                             | mismatch for V gene.               |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-gap-open-penalty``    | Penalty (> 0) to open a gap for    | ``6``          |
|                             | V gene.                            |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-gap-extend-penalty``  | Penalty (> 0) to extend a gap for  | ``6``          |
|                             | V gene.                            |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-wordsize``            | Word size to find hotspots by      | ``10``         |
|                             | BLAST for V gene.                  |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-evalue-cutoff``       | E-value                            | ``1e-50``      |
|                             | threshold for assigning V gene.    |                |   
+-----------------------------+------------------------------------+----------------+
| ``--d-blastdb``             | Path to BLAST database that        |                |
|                             | created from the FASTA file of     |                |
|                             | ``-d``.                            |                | 
+-----------------------------+------------------------------------+----------------+
| ``--d-match-score``         | Score (> 0) for a nucleotide match | ``1``          |
|                             | for D gene.                        |                |
+-----------------------------+------------------------------------+----------------+
| ``--d-mismatch-score``      | Score (< 0) for a nucleotide       | ``-1``         |
|                             | mismatch for D gene.               |                |
+-----------------------------+------------------------------------+----------------+
| ``--d-gap-open-penalty``    | Penalty (> 0) to open a gap for    | ``0``          |
|                             | D gene.                            |                |
+-----------------------------+------------------------------------+----------------+
| ``--d-gap-extend-penalty``  | Penalty (> 0) to extend a gap for  | ``2``          |
|                             | D gene.                            |                |
+-----------------------------+------------------------------------+----------------+
| ``--d-wordsize``            | Word size to find hotspots by      | ``4``          |
|                             | BLAST for D gene.                  |                |
+-----------------------------+------------------------------------+----------------+
| ``--d-evalue-cutoff``       | E-value                            | ``1e-2``       |
|                             | threshold for assigning D gene.    |                |   
+-----------------------------+------------------------------------+----------------+
| ``--j-blastdb``             | Path to BLAST database that        |                |
|                             | created from the FASTA file of     |                |
|                             | ``-j``.                            |                | 
+-----------------------------+------------------------------------+----------------+
| ``--j-match-score``         | Score (> 0) for a nucleotide match | ``3``          |
|                             | for J gene.                        |                |
+-----------------------------+------------------------------------+----------------+
| ``--j-mismatch-score``      | Score (< 0) for a nucleotide       | ``-3``         |
|                             | mismatch for J gene.               |                |
+-----------------------------+------------------------------------+----------------+
| ``--j-gap-open-penalty``    | Penalty (> 0) to open a gap for    | ``6``          |
|                             | J gene.                            |                |
+-----------------------------+------------------------------------+----------------+
| ``--j-gap-extend-penalty``  | Penalty (> 0) to extend a gap for  | ``6``          |
|                             | J gene.                            |                |
+-----------------------------+------------------------------------+----------------+
| ``--j-wordsize``            | Word size to find hotspots by      | ``7``          |
|                             | BLAST for J gene.                  |                |
+-----------------------------+------------------------------------+----------------+
| ``--j-evalue-cutoff``       | E-value                            | ``1e-5``       |
|                             | threshold for assigning J gene.    |                |   
+-----------------------------+------------------------------------+----------------+
| ``--v-motif``               | The regular expression pattern     | ``YYC``        |
|                             | of the motif on V gene to identify |                |
|                             | CDR3 segment.                      |                |
|                             | (e.g, ``YYC`` ``[FY]YC`` )         |                |   
+-----------------------------+------------------------------------+----------------+
| ``--j-motif``               | The regular expression pattern     | ``WG.G``       |
|                             | of the motif on J gene to identify |                |
|                             | CDR3 segment.                      |                |
|                             | (e.g, ``WG.G`` ``FG.G`` ``G.G``)   |                |   
+-----------------------------+------------------------------------+----------------+



Output
^^^^^^


The :command:`parse` action generates some intermediate results with TSV format
and the final results with PYDAIR format.



+--------------------------+------------+-----------------------------------------------------+
| File name                | Format     | File contents                                       |
+==========================+============+=====================================================+
| *<prefix>*.v.blast.txt   | TSV        | BLAST results for V gene.                           |
+--------------------------+------------+-----------------------------------------------------+
| *<prefix>*.j.blast.txt   | TSV        | BLAST results for J gene.                           |
+--------------------------+------------+-----------------------------------------------------+
| *<prefix>*.unaligned.fa  | FASTA      | Unaligned region sequences that was used for        |
|                          |            | BLAST to identify D.                                |
+--------------------------+------------+-----------------------------------------------------+
| *<prefix>*.d.blast.txt   | TSV        | BLAST results for D gene.                           |
+--------------------------+------------+-----------------------------------------------------+
| *<prefix>*.vdj.pydair    | PYDAIR     | The final results of :command:`parse` mode.         |
+--------------------------+------------+-----------------------------------------------------+




stats action mode
=================

The :command:`stats` action mode summarizes the identification results
and creates an HTML report.


Usage
^^^^^

.. code-block:: bash
    
    pydair stats -i output1.pydair output2.pydair output3.pydair  \
                 -n Sample1 Sample2 Sample3                       \
                 -o stats_result

The three options, i.e., ``-i``, ``-n``, and ``-o``
are required options.


Option
^^^^^^

+--------------------------------+------------------------------------+----------------+
| Options                        | Descriptions                       | Default        |
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
| ``--discard-ambiguous-D``      | If ``False``, summarize all        | ``False``      |
|                                | sequences regardless the D segment |                |
|                                | is identified or not.              |                |
|                                | If ``True``, the summarization is  |                |
|                                | performed after discarding         |                |
|                                | sequences with ambiguous D         |                |
|                                | segment.                           |                |
+--------------------------------+------------------------------------+----------------+
| ``--productive-only``          | If ``False``, summarize all        | ``False``      |
|                                | sequences regardless productive or |                |
|                                | nonproductive sequences.           |                |
|                                | If ``True``, summarize only the    |                |
|                                | productive sequences.              |                |
+--------------------------------+------------------------------------+----------------+
| ``--estimate-vdj-combination`` | If ``True``, perform rarefaction   | ``False``      |
|                                | analysis to study the diversity of |                |
|                                | VDJ combinations.                  |                |
+--------------------------------+------------------------------------+----------------+




Output
^^^^^^


+-------------------------------------+------------+-----------------------------------------------------+
| File name                           | Format     | File contents                                       |
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
    QN Seq-0-25045;
    VN J03617|IGHV3-53*03|Homo
    DN X13972|IGHD5-12*01|Homo
    JN J00256|IGHJ3*01|Homo
    OP .
    QA           CCGTGGAGTCTGGAGGAGGCTTGATCCAGCCTGAGGGGTCCCTGAGACTCTCCCATGCAGCCTCTGGGTTCACTGTCAGTAGAAACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGTCTTCTATTTATAGCGGTGGTAGCACATACTACGCAGACTCTGTGAAGGGCCGATTCACCATCTCCTGAGACTATTCCAAGAACACGCTGTATCTTCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCTAGAACTATAGTGGCTACGATTTTTTTATGACTGGGGCCAAGGGACAATGGTCAC
    VA GAGGTGCAGCTGGTGGAGTCTGGAGGAGGCTTGATCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGGTTCACCGTCAGTAGCAACTACATGAGCTGGGTCCGCCAGCCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGT----TATTTATAGCGGTGGTAGCACATACTACGCAGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTTCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCTAGGGA
    JA                                                                                                                                                                                                                                                                                                                    TGATGCTTTTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
    UA                                                                                                                                                                                                                                                                                                        ACTATAGTGGCTACGATT
    CA                                                                                                                                                                                                                                                                                                  GCTAGAACTATAGTGGCTACGATTTTTTTATGAC
    VD GGA
    JD TGATGC
    VJ AACTATAGTGGCTACGATT
    #CDR3AA ARTIVATIFL*
    #AL QSTART  QEND    SSTART  SEND    IDENTITY    SCORE
    AL QV 3 284 13  290 96.099  7.10e-112
    AL QD . .   .   .   100.000 0.003
    AL QJ 304   336 7   39  93.939  1.34e-09
    AL QU 285   303 .   .   .   .
    AL QC 279   313 .   .   .   .
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
|           | ``.`` indicates that the sequence is           |
|           | unproductive.                                  |
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
| ``VD``    | Deleted bases at the 5'-end of V gene.         |
+-----------+------------------------------------------------+
| ``JD``    | Deleted bases at the 3'-end of J gene.         |
+-----------+------------------------------------------------+
| ``VJ``    | Inserted bases between 5'-end of V and 3'-end  |
|           | of J gene. This segment contains D segment,    |
|           | and this segment is containted in CDR3 segment.|
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




sim action mode
===============

The :command:`sim` action mode generates artificial IgH sequences
and outputs as a FASTA file.

Usage
^^^^^

.. code-block:: bash
    
    pydair sim -n 10000 -o output_sim.fa \
               --v-fasta v.fa \
               --d-fasta d.fa \
               --j-fasta j.fa
    

The option `-n` is used for specifying the number of sequences should be generated,
`-o` is used for specifying file path to save the generated sequences.
In addition, `--v-fasta`, `--d-fasta`, and `--j-fasta` is used for
specifiying the population pools for VDJ ramdom sampling.


Options
^^^^^^^

In addition to the required options,
the parameters of nucleotide additions, deletions and mutations can be specified.

+-----------------------------+------------------------------------+----------------+
| Options                     | Description                        | Default        |
+=============================+====================================+================+
| ``-n``                      | The number of sequences should be  | `10000`        |
|                             | generated.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``-o``                      | Path to file for writting results. |                |
+-----------------------------+------------------------------------+----------------+
| ``--v-fasta``               | Path to FASTA file of V gene       |                |
|                             | sequence (the population for       |                |
|                             | sampling).                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--d-fasta``               | Path to FASTA file of D gene       |                |
|                             | sequence (the population for       |                |
|                             | sampling).                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--j-fasta``               | Path to FASTA file of J gene       |                |
|                             | sequence (the population for       |                |
|                             | sampling).                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-v-5del``              | The mean of Poisson distribution   | ``10``         |
|                             | that for sampling the number of    |                |
|                             | nucleotides deleting from 5'-end   |                |
|                             | of V gene.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-v-3del``              | The mean of Poisson distribution   | ``3``          |
|                             | that for sampling the number of    |                |
|                             | nucleotides deleting from 5'-end   |                |
|                             | of V gene.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-d-5del``              | The mean of Poisson distribution   | ``3``          |
|                             | that for sampling the number of    |                |
|                             | nucleotides deleting from 5'-end   |                |
|                             | of D gene.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-d-3del``              | The mean of Poisson distribution   | ``3``          |
|                             | that for sampling the number of    |                |
|                             | nucleotides deleting from 5'-end   |                |
|                             | of D gene.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-j-5del``              | The mean of Poisson distribution   | ``5``          |
|                             | that for sampling the number of    |                |
|                             | nucleotides deleting from 5'-end   |                |
|                             | of J gene.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-j-3del``              | The mean of Poisson distribution   | ``10``         |
|                             | that for sampling the number of    |                |
|                             | nucleotides deleting from 5'-end   |                |
|                             | of J gene.                         |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-vd-ins``              | The mean of Poisson distribution   | ``5``          |
|                             | that for sampling the number of    |                |
|                             | nucleotides inserting into VD      |                |
|                             | junciton.                          |                |
+-----------------------------+------------------------------------+----------------+
| ``--n-dj-ins``              | The mean of Poisson distribution   | ``5``          |
|                             | that for sampling the number of    |                |
|                             | nucleotides inserting into DJ      |                |
|                             | junciton.                          |                |
+-----------------------------+------------------------------------+----------------+
| ``--p-mutation``            | The probability to mutate a        | ``0.05``       |
|                             | nucleotide.                        |                |
+-----------------------------+------------------------------------+----------------+






eval action mode
================

The :command:`eval` action mode is used for evaluating performances of PyDAIR.
It usually used after ``sim`` and ``parse`` mode.

.. code-block:: bash
    
    pydair sim -n 10000 -o output_sim.fa \
               --v-fasta v.fa \
               --d-fasta d.fa \
               --j-fasta j.fa
    
    pydair parse -q output_sim.fa \
                 -o output1                \
                 -v v.fa                   \
                 -d d.fa                   \
                 -j j.fa                   \
                 --v-blastdb blastdb_v     \
                 --d-blastdb blastdb_d     \
                 --j-blastdb blastdb_j 
    
    pydair eval -o eval_result.txt
                --sim-condition output_sim.fa \
                --parse-result output1.vdj.pydair
    

The results are saved into text file with TSV format.


Options
^^^^^^^

+-----------------------------+------------------------------------+----------------+
| Options                     | Description                        | Default        |
+=============================+====================================+================+
| ``-o``                      | Path to file for writting results. |                |
+-----------------------------+------------------------------------+----------------+
| ``--sim-condition``         | The FASTA file which generated by  |                |
|                             | `sim` action mode.                 |                |
+-----------------------------+------------------------------------+----------------+
| ``--parse-result``          | The PYDAIR file which generated by |                |
|                             | `parse` action mode.               |                |
+-----------------------------+------------------------------------+----------------+



