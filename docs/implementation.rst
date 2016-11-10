==============
Implementation
==============


PyDAIR offers two commands, :command:`pydair parse` and :command:`pydair stats`,
to perform repertoire diversity study.
The :command:`pydair parse` command is used for identification of V, D, and J genes
and determination of CDR3 sequences from each of Ig-Seq sequence,
and outputting results into a flat file.
The :command:`pydair stats` command is used for summarizing results and
outputting into tab-delimited text files.



Rep-Seq sequence parse
======================

:command:`pydair parse` comamnd is used to identify V, D, and J genes
that usaged in query sequences, and determine CDR3 in sequence.
It can be explained in fifth steps.

#. To identify V gene that used in Rep-Seq IgH sequence,
   an IgH sequences is aligned to V gene database using BLAST.
   The best hit of V gene in database is assigned to V gene that
   used in IgH sequence.
   Parameters of BLAST can be changed by users.
#. Same to assignment of V gene, the J gene is identified by
   aligning IgH sequence to J gene database using BLAST.
#. To identify CDR3 sequence, regular expression matching is
   used to seek the YYC motif from the end of V-aligned region
   and seek the WGxG motif from the start of J-aligned region.
   This is because CDR3 sequence is defined between these two
   highly conserved motifs.
#. To identify D gene, the identified CDR3 sequence is aligned
   to D gene database using BLAST.
   However, if CDR3 sequence is unidentifiable,
   the region that impossible to align to V and J genes
   (unaligned region) is used instead of CDR3 sequence.
#. To identify the open reading frame (ORF) of IgH sequence,
   regular expression matching is used to seek the region
   beginning with a methionine in V region and without any
   stop codon in downstream, since the ORF is defined by the
   start codon in V region of IgH.


.. image:: images/pydair-parse.png
    :align: center
    




Repertoire diversity study
==========================


The :command:`pydair stats` command is to calculate the frequencies of
V, D, and J gene usages, and the distribution of CDR3 sequence
length with tab-delimited text file.




PYDAIR format
=============

PYDAIR format is an human readable text file format.
Each entry in PYDAIR format file represents a result of
VDJ identification and CDR3 sequence determination for a query sequence. 
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

+------------+-------------------------------------------------+
| attribute  | definition                                      |
+============+=================================================+
| ``QV``     | Alignment and BLAST results between query and V |
|            | gene.                                           |
+------------+-------------------------------------------------+
| ``QD``     | Alignment and BLAST results between query and V |
|            | gene.                                           |
+------------+-------------------------------------------------+
| ``QJ``     | Alignment and BLAST results between query and V |
|            | gene.                                           |
+------------+-------------------------------------------------+
| ``QU``     | Alignment resulTs between query and unaligned   |
|            | sequence.                                       |
+------------+-------------------------------------------------+
| ``QC``     | Alignment resulTs between query and CDR3        |
|            | sequence.                                       |
+------------+-------------------------------------------------+

Each attribute contains the five columns that separated by TAB.
The five columns indicates that the alignement start and end
position in query sequence,
the alignemnt start and end position in subject (V, D, J, unaligned, and CDR3) sequence,
and the identity and BLAST score of alignemnt.





