*Python library for diversity analysis for immune repertoire*

======
PyDAIR
======



Immunoglobulin (Ig) or antibody repertoire is one of the key players in
a potent adaptive immune system in animals.
It is estimated that Ig has trillion order of different sequences in human,
which enable to recognize various pathogens and give rise to specific immune
responses continuously.

The diversity of Ig is mainly derived from VDJ recombination on Ig heavy chains (IgH).
The recombination event process includes:
(i) single of V (variable) gene, D (diversity) gene,
and J (joining) gene are randomly selected from sets of germline gene to rearrangement;
(ii) random numbers of nucleotides are trimmed from both end of V, D and J genes;
and (iii) random numbers of nucleotides are inserted into V-D and D-J junctions.
The region between the V-D and D-J junctions on IgH,
the complementarity-determining region 3 (CDR3),
is the most variable component of the IgH and is a predominant determinant of binding specificity.

Analyzing of diversity of immune repertoire gives us an opportunity to understand adaptive immunity.
Deep profiling of Ig nucleotide sequences by means of high-throughput sequencing,
termed as Ig-Seq, has become an attractive way for studying Ig repertoire diversity at high resolution.
In general, Ig-Seq outputs millions of “reads”,
and each read represents single of IgH (or IgL) nucleotide sequences.
Thus, the beginning of step to study Ig repertoire diversity involves
identification of V, D, and J segments for each read.

PyDAIR is a tool that aims to study diversity of IgH sequences based on Ig-Seq data.
PyDAIR identifies V, D, and J gene segments using BLAST,
and identifies CDR3 segment using regular expression matching from millions of reads.
Furthermore, PyDAIR enables detection of insertions and deletions on VD and DJ junctions,
identification of CDR3 segment, diversity measurement of VDJ combinations,
and results visualization in HTML report.
PyDAIR can be applied to any organism without database restriction.





.. toctree::
    :hidden:
    :numbered:
    
    immunology
    implementation
    usage
    casestudies
    
    
