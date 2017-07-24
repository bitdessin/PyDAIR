======
Usages
======


Installation
============


PyDAIR is a Python package which can be installed via :command:`pip` command.

.. code-block:: text

    pip install pydair


In addtion to installation of PyDAIR,
NCBI BLAST+ is also required to be installed.
For Linux users, BALST can be installed via :command:`apt-get` command.


.. code-block:: text
    
    apt-get install ncbi-blast+


For Macintosh users, BLAST can be installed via :command:`brew` command
after adding :command:`homebrew/science` repository.


.. code-block:: text
    
    brew tap homebrew/science
    brew install blast







parse mode
==========

The :command:`parse` mode assigns V, D, J, and CDR3 segments,
and detects the deletions and insertions in IgH sequences.

Usage
^^^^^

.. code-block:: text
    
    pydair parse -q input_igh_sequences.fa \
                 -o output1                \
                 -v v.fa                   \
                 -d d.fa                   \
                 -j j.fa                   \
                 --v-blastdb blastdb_v     \
                 --d-blastdb blastdb_d     \
                 --j-blastdb blastdb_j 


The eight options, i.e., :command:`-q`, :command:`-o`, :command:`-v`,
:command:`-d`, :command:`-j`, :command:`--v-blastdb`, :command:`--d-blastdb`,
and :command:`--j-blastdb`, are required options.


Options
^^^^^^^

In addition to the required options, 
the other options such as BLAST parameters for VDJ identification can be specified.
All options used in :command:`parse` mode are shown in the following table.


+------------------------------------+------------------------------------+----------------+
| Options                            | Description                        | Default        |
+====================================+====================================+================+
| :command:`-q`                      | Path to FASTA file that contains   |                |
|                                    | IgH seqeunces.                     |                |
+------------------------------------+------------------------------------+----------------+
| :command:`-o`                      | Path to file for writting results. |                |
+------------------------------------+------------------------------------+----------------+
| :command:`-v`                      | Path to FASTA file of V gene       |                |
|                                    | sequence.                          |                |
+------------------------------------+------------------------------------+----------------+
| :command:`-d`                      | Path to FASTA file of D gene       |                |
|                                    | sequence.                          |                |
+------------------------------------+------------------------------------+----------------+
| :command:`-j`                      | Path to FASTA file of J gene       |                |
|                                    | sequence.                          |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-blastdb`             | Path to BLAST database that        |                |
|                                    | created from the FASTA file of     |                |
|                                    | :command`-v`.                      |                | 
+------------------------------------+------------------------------------+----------------+
| :command:`--v-match-score`         | Score (> 0) for a nucleotide match | ``3``          |
|                                    | for V gene.                        |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-mismatch-score`      | Score (< 0) for a nucleotide       | ``-3``         |
|                                    | mismatch for V gene.               |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-gap-open-penalty`    | Penalty (> 0) to open a gap for    | ``6``          |
|                                    | V gene.                            |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-gap-extend-penalty`  | Penalty (> 0) to extend a gap for  | ``6``          |
|                                    | V gene.                            |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-wordsize`            | Word size to find hotspots by      | ``10``         |
|                                    | BLAST for V gene.                  |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-evalue-cutoff`       | E-value                            | ``1e-50``      |
|                                    | threshold for assigning V gene.    |                |   
+------------------------------------+------------------------------------+----------------+
| :command:`--d-blastdb`             | Path to BLAST database that        |                |
|                                    | created from the FASTA file of     |                |
|                                    | :command:`-d`.                     |                | 
+------------------------------------+------------------------------------+----------------+
| :command:`--d-match-score`         | Score (> 0) for a nucleotide match | ``1``          |
|                                    | for D gene.                        |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--d-mismatch-score`      | Score (< 0) for a nucleotide       | ``-1``         |
|                                    | mismatch for D gene.               |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--d-gap-open-penalty`    | Penalty (> 0) to open a gap for    | ``0``          |
|                                    | D gene.                            |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--d-gap-extend-penalty`  | Penalty (> 0) to extend a gap for  | ``2``          |
|                                    | D gene.                            |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--d-wordsize`            | Word size to find hotspots by      | ``4``          |
|                                    | BLAST for D gene.                  |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--d-evalue-cutoff`       | E-value                            | ``1e-2``       |
|                                    | threshold for assigning D gene.    |                |   
+------------------------------------+------------------------------------+----------------+
| :command:`--j-blastdb`             | Path to BLAST database that        |                |
|                                    | created from the FASTA file of     |                |
|                                    | :command:`-j`.                     |                | 
+------------------------------------+------------------------------------+----------------+
| :command:`--j-match-score`         | Score (> 0) for a nucleotide match | ``3``          |
|                                    | for J gene.                        |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--j-mismatch-score`      | Score (< 0) for a nucleotide       | ``-3``         |
|                                    | mismatch for J gene.               |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--j-gap-open-penalty`    | Penalty (> 0) to open a gap for    | ``6``          |
|                                    | J gene.                            |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--j-gap-extend-penalty`  | Penalty (> 0) to extend a gap for  | ``6``          |
|                                    | J gene.                            |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--j-wordsize`            | Word size to find hotspots by      | ``7``          |
|                                    | BLAST for J gene.                  |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--j-evalue-cutoff`       | E-value                            | ``1e-5``       |
|                                    | threshold for assigning J gene.    |                |   
+------------------------------------+------------------------------------+----------------+
| :command:`--v-motif`               | The regular expression pattern     | ``YYC``        |
|                                    | of the motif on V gene to identify |                |
|                                    | CDR3 segment.                      |                |
|                                    | (e.g, ``YYC`` and ``[FY]YC`` )     |                |   
+------------------------------------+------------------------------------+----------------+
| :command:`--j-motif`               | The regular expression pattern     | ``WG.G``       |
|                                    | of the motif on J gene to identify |                |
|                                    | CDR3 segment.                      |                |
|                                    | (e.g, ``WG.G``, ``FG.G``, and      |                |
|                                    | ``[WF]G.G``)                       |                |   
+------------------------------------+------------------------------------+----------------+



Output
^^^^^^


The :command:`parse` mode generates some intermediate results with TSV format
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
| *<prefix>*.vj.pydair     | PYDAIR     | Intermediate result which lacks of results of D     |
|                          |            | gene identification.                                |
+--------------------------+------------+-----------------------------------------------------+
| *<prefix>*.vdj.pydair    | PYDAIR     | The final results of :command:`parse` mode.         |
+--------------------------+------------+-----------------------------------------------------+

Except for the :file:`<prefix>.vdj.pydair` file, other files can be removed to save disk space.




stats mode
==========

The :command:`stats` mode summarizes the identification results and creates an HTML report.


Usage
^^^^^

.. code-block:: text
    
    pydair stats -i output1.pydair output2.pydair output3.pydair  \
                 -n Sample1 Sample2 Sample3                       \
                 -o stats_result


The three options, i.e., :command:`-i`, :command:`-n`,
and :command:`-o`, are required options.


Option
^^^^^^

+---------------------------------------+------------------------------------+----------------+
| Options                               | Descriptions                       | Default        |
+=======================================+====================================+================+
| :command:`-i`                         | Path to ``PYDAIR`` format files.   |                |
|                                       | Multiple files should be separated |                |
|                                       | by a blank.                        |                |
+---------------------------------------+------------------------------------+----------------+
| :command:`-n`                         | Sample names of each ``PYDAIR``    |                |
|                                       | files. Multiple names should be    |                |
|                                       | separated by a blank.              |                |
+---------------------------------------+------------------------------------+----------------+
| :command:`-o`                         | Prefix for writting results.       |                |
+---------------------------------------+------------------------------------+----------------+
| :command:`--discard-ambiguous-D`      | If ``False``, summarize all        | ``False``      |
|                                       | sequences regardless the D segment |                |
|                                       | is identified or not.              |                |
|                                       | If ``True``, the summarization is  |                |
|                                       | performed after discarding         |                |
|                                       | sequences with ambiguous D         |                |
|                                       | segment.                           |                |
+---------------------------------------+------------------------------------+----------------+
| :command:`--productive-only`          | If ``False``, summarize all        | ``False``      |
|                                       | sequences regardless productive or |                |
|                                       | nonproductive sequences.           |                |
|                                       | If ``True``, summarize only the    |                |
|                                       | productive sequences.              |                |
+---------------------------------------+------------------------------------+----------------+
| :command:`--estimate-vdj-combination` | If ``True``, perform rarefaction   | ``False``      |
|                                       | analysis to study the diversity of |                |
|                                       | VDJ combinations.                  |                |
+---------------------------------------+------------------------------------+----------------+




Output
^^^^^^

The summarization results are saved in TSV format files and an HTML file.
Users can open the HTML file with web browsers to check the summarization.

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







sim mode
===============

The :command:`sim` mode generates artificial IgH sequences and outputs as a FASTA file.

Usage
^^^^^

.. code-block:: text
    
    pydair sim -n 10000 -o output_sim.fa \
               --v-fasta v.fa \
               --d-fasta d.fa \
               --j-fasta j.fa
    

The option :command:`-n` is used for specifying the number of sequences should be generated,
:command:`-o` is used for specifying file path to save the generated sequences.
In addition, :command:`--v-fasta`, :command:`--d-fasta`, and :command:`--j-fasta` are used for
specifiying the population pools for VDJ ramdom sampling.


Options
^^^^^^^

In addition to the required options,
the parameters of nucleotide additions, deletions and mutations can be specified.

+------------------------------------+------------------------------------+----------------+
| Options                            | Description                        | Default        |
+====================================+====================================+================+
| :command:`-n`                      | The number of sequences should be  | ``10000``      |
|                                    | generated.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`-o`                      | Path to file for writting results. |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--v-fasta`               | Path to FASTA file of V gene       |                |
|                                    | sequence (the population for       |                |
|                                    | sampling).                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--d-fasta`               | Path to FASTA file of D gene       |                |
|                                    | sequence (the population for       |                |
|                                    | sampling).                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--j-fasta`               | Path to FASTA file of J gene       |                |
|                                    | sequence (the population for       |                |
|                                    | sampling).                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-v-5del`              | The mean of Poisson distribution   | ``10``         |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides deleting from 5'-end   |                |
|                                    | of V gene.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-v-3del`              | The mean of Poisson distribution   | ``3``          |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides deleting from 3'-end   |                |
|                                    | of V gene.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-d-5del`              | The mean of Poisson distribution   | ``3``          |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides deleting from 5'-end   |                |
|                                    | of D gene.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-d-3del`              | The mean of Poisson distribution   | ``3``          |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides deleting from 3'-end   |                |
|                                    | of D gene.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-j-5del`              | The mean of Poisson distribution   | ``5``          |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides deleting from 5'-end   |                |
|                                    | of J gene.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-j-3del`              | The mean of Poisson distribution   | ``10``         |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides deleting from 3'-end   |                |
|                                    | of J gene.                         |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-vd-ins`              | The mean of Poisson distribution   | ``5``          |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides inserting into VD      |                |
|                                    | junciton.                          |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--n-dj-ins`              | The mean of Poisson distribution   | ``5``          |
|                                    | that for sampling the number of    |                |
|                                    | nucleotides inserting into DJ      |                |
|                                    | junciton.                          |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--p-mutation`            | The probability to mutate a        | ``0.05``       |
|                                    | nucleotide.                        |                |
+------------------------------------+------------------------------------+----------------+






eval mode
================

The :command:`eval` mode is used for evaluating performances of PyDAIR.
This mode usually is used after :command:`sim` and :command:`parse` modes.


.. code-block:: text
    
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
    

The evaluation results are saved into text file with TSV format.


Options
^^^^^^^

+------------------------------------+------------------------------------+----------------+
| Options                            | Description                        | Default        |
+====================================+====================================+================+
| :command:`-o`                      | Path to file for writting results. |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--sim-condition`         | The FASTA file which generated by  |                |
|                                    | :command:`sim` mode.               |                |
+------------------------------------+------------------------------------+----------------+
| :command:`--parse-result`          | The PYDAIR file which generated by |                |
|                                    | :command:`parse` mode.             |                |
+------------------------------------+------------------------------------+----------------+






PYDAIR format
=============

PYDAIR format is an human readable text file format.
PYDAIR format file can contain multiple entries.
Each entry in PYDAIR format file represents a result of
V, D, J and CDR3 segment identification of a query IgH sequence. 
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
The line begining with a ``#`` is a comment line,
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

