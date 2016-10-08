======
Usages
======

Installation
============


The simplest way to install the latest PyDAIR is to use :command:`pip` command.

.. code-block:: bash

    pip install --user

PyDAIR requires following software to be installed.

* Python 2.7 or Python 3.4
* `NumPy <http://www.numpy.org/>`_
* `Pandas <http://pandas.pydata.org/>`_
* `matplotlib <http://matplotlib.org/>`_
* `BioPython <http://biopython.org/>`_
* `NCBI BLAST+ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_



Command and options
===================

PyDAIR command (:command:`pydair`) consists of the two actions, :command:`parse` and :command:`stats`.


+-------------------+--------------------------------------------------------------+
| actions           | function                                                     |
+===================+==============================================================+
| :command:`parse`  | Parse IgH sequences. Namely, identification of V, D, and J   |
|                   | genes, and determination of CDR3 sequences.                  |
+-------------------+--------------------------------------------------------------+
| :command:`stats`  | Summarize the parsed results and save them into TSV file.    |
+-------------------+--------------------------------------------------------------+



parse action
^^^^^^^^^^^^

:command:`parse` action is to identify V, D, and J genes from each IgH sequence
by aligning IgH sequence to germline (V, D, and J) database using NCBI BLAST+.
It requires IgH sequences (in FASTA format),
germline sequences (in FASTA format),
BLAST databases of germiline sequences,
and BLAST parameters.

PyDAIR generates several files to save the intermediate results,
such as BLAST results, region that cannot be aligned to V and J genes.
The final result is saved into :file:`output1.pydair` file.
If there several samples, :command:`parse` action should be run several times for each sample.


.. code-block:: bash
    
    pydair parse -q input_igh_sequences.fa \
                 -v v.fa                   \
                 -d d.fa                   \
                 -j j.fa                   \
                 --v-blastdb blastdb_v     \
                 --d-blastdb blastdb_d     \
                 --j-blastdb blastdb_j     \
                 -o output1


+----------------------------+------------------------------------+----------------+
| options                    |                                    | defult         |
+============================+====================================+================+
| ``-q``                     | Path to FASTA file that contains   |                |
|                            | IgH seqeunces.                     |                |
+----------------------------+------------------------------------+----------------+
| ``-o``                     | Path to file for writting results. |                |
+----------------------------+------------------------------------+----------------+
| ``-f``                     | File format of results. ``pydair`` | ``pydair``     |
|                            | or ``simple`` can be specified.    |                |
+----------------------------+------------------------------------+----------------+
| ``-s``                     | Species.                           | ``fugu``       |
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
|                            | created from FASTA file of  ``v``. |                | 
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
| ``-v-evalue-cutoff``       | Expectation value (e-value)        | ``1e-10``      |
|                            | threshold for assiging V gene.     |                |   
+----------------------------+------------------------------------+----------------+
| ``-d-blastdb``             | Path to BLAST database that        |                |
|                            | created from FASTA file of  ``d``. |                | 
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
| ``-d-evalue-cutoff``       | Expectation value (e-value)        | ``1e-2``       |
|                            | threshold for assiging D gene.     |                |   
+----------------------------+------------------------------------+----------------+
| ``-j-blastdb``             | Path to BLAST database that        |                |
|                            | created from FASTA file of ``j``.  |                | 
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
| ``-j-evalue-cutoff``       | Expectation value (e-value)        | ``1e-5``       |
|                            | threshold for assiging J gene.     |                |   
+----------------------------+------------------------------------+----------------+







stats action
^^^^^^^^^^^^

The statistical summaries are calculated by :command:`stats` action.


.. code-block:: bash
    
    pydair stats -i output1.pydair output2.pydair output3.pydair  \
                 -n Fugu1 Fugu2 Fugu3                             \
                 -o stats_result                                  \
                 --contain_ambiguous_D



+--------------------------------+------------------------------------+----------------+
| options                        |                                    | defult         |
+================================+====================================+================+
| ``-i``                         | Path to ``pydair`` format files.   |                |
|                                | Multiple files should be separated |                |
|                                | by a blank.                        |                |
+--------------------------------+------------------------------------+----------------+
| ``-n``                         | Sample names of each ``pydair``    |                |
|                                | files. Multiple names should be    |                |
|                                | separated by a blank.              |                |
+--------------------------------+------------------------------------+----------------+
| ``-o``                         | Prefix for writting results.       |                |
+--------------------------------+------------------------------------+----------------+
| ``--contain_ambiguous_D``      | Contain IgH sequence that has      | ``True``       |
|                                | unidentifiable D genes.            |                |
+--------------------------------+------------------------------------+----------------+
| ``--containe_stopcodon``       | Contain IgH sequence that has top  | ``False``      |
|                                | codons.                            |                |
+--------------------------------+------------------------------------+----------------+
| ``--figure-format``            | Figure format.                     | ``png``        |
+--------------------------------+------------------------------------+----------------+
| ``--figure-dpi``               | Figure DPI.                        | ``300``        |
+--------------------------------+------------------------------------+----------------+
| ``--figure-style``             | Figure style. 'classic', 'ggplot', | ``ggplto``     |
|                                | can be specified.                  |                |
+--------------------------------+------------------------------------+----------------+
| ``--estimate-vdj-combination`` | Rarefaction study for              | ``False``      |
|                                | VDJ combinations.                  |                |
+--------------------------------+------------------------------------+----------------+



