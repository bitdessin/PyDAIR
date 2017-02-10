==============
Implementation
==============


PyDAIR offers two action modes, :command:`parse` and :command:`stats`.
The :command:`parse` action is used to identify V, D, J and CDR3 segments in IgH sequences with given Ig-Seq data.
On the other hand, the :command:`stats` action is used to summarize and visualize the analysis results.



parse action
============

:command:`parse` command is used to identify V, D, J, and CDR3 segments,
and detect the indels around V-D and D-J junctions.

#. **V segment identification** To identify V segment,
   PyDAIR aligns each IgH sequence against user-defined BLAST database of V gene with :command:`blastn`.
   The best matching V gene is selected.
#. **J segment identification** To identify J segment,
   PyDAIR aligns each IgH sequence against user-defined BLAST database of J gene with :command:`blastn`.
   The best matching J gene is selected.
#. **CDR3 segment identification** To identify CDR3 segment,
   PyDAIR search for YYC and WGxG motifs using regular expression matching.
#. **Indels detection**  To identify insertion and deletion on V-D and D-J junctions,
   PyDAIR calculates the number of unaligned bases from the end of V and the start of J genes.
#. **D segment identification** To identify D segment,
   PyDAIR aligns the each insertion segment sequence against user-defined BLAST database of D gene with :command:`blastn`.
   The best matching D gene is selected.
#. **Open reading frame detection**
   PyDAIR searches for the start codon (ATG) and the stop codons (TAG, TAA, and TGA).
   PyDAIR tries the three possible reading frames to determine ORF.
   The sequence is defined as *productive* if one of reading frames dose not contain any stop codons.
   Otherwise, the sequence is defined as *unproductive*.


.. image:: images/pydair-parse.png
    :align: center
    


stats action
============

:command:`stats` command is used for summarizing V, D, and J segment usage frequencies,
and the distribution of length of CDR3 segment.
The summarized results are saved into TSV (tab-delimited) files.
In addition, the summarized results are visualized with some charts in an HTML report.



