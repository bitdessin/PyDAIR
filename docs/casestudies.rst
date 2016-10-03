============
Case studies
============

.. note:: Make sure PyDIAR and BLAST have been already installed.




Simulation study with Joinsim dataset
=====================================

Artificial IgH sequences, human germline gene databases,
and scripts that used in the simulation study are available on
`GitHub PyDAIR repository <https://github.com/bioinfoteam/PyDAIR>`_.
We use ``git clone`` command to download these data from GitHub.


.. code-block:: bash
    
    git clone git@github.com:bioinfoteam/PyDAIR.git


The data are saved in `PyDAIR/casestudies/joinsim`.
We use ``cd`` command to go to `joinsim` directory.


.. code-block:: bash
    
    cd PyDAIR/casestudies/joinsim




Preparation
^^^^^^^^^^^

Before analysis, we create BLAST database with human
germline gene sequences using ``makeblastdb``.


.. code-block:: bash
    
    cd db
    makeblastdb -in human.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


The artificial sequences obtained from JoinSimulation,
is saved in `data/joinsim.txt` directory with TSV format.
We first extract sequences with ORF from whole data set
to create a subset,
and save them to `data/joinsim.sub.txt`.
Then, we convert the whole data set and subset into FASTA format file.


.. code-block:: bash
    
    grep -v "out of frame" data/joinsim.txt | grep -v "stop codon"  > data/joinsim.sub.txt
    python ./bin/convert_csv_to_fastq.py ./data/joinsim.txt ./data/joinsim.fa
    python ./bin/convert_csv_to_fastq.py ./data/joinsim.sub.txt ./data/joinsim.sub.fa






Analysis
^^^^^^^^

We use ``pydair parse`` command to assign VDJ genes,
and determine CDR3 sequences for whole data set and subset, respectively.


.. code-block:: bash
    
    pydair parse -s human -q data/joinsim.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o results/joinsim
    
    pydair parse -s human -q data/joinsim.sub.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o results/joinsim.sub


Then, we use the original Python scripts to calculate
the number of sequences that are correctly and incorrectly asssigned.


.. code-block:: bash
    
    python ./bin/calc_accuracy.py ./data/joinsim.txt \
                                  ./results/joinsim.vdj.pydair.simple \
                                  ./results/joinsim.stats.txt
    
    python ./bin/calc_accuracy.py ./data/joinsim.sub.txt \
                                  ./results/joinsim.sub.vdj.pydair.simple \
                                  ./results/joinsim.sub.stats.txt


The calculation results are saved into 
`joinsim.stats.txt` for whole data set,
and `joinsim.sub.stats.txt` for subset.

Finally, we use ``pydair stats`` commands to create
the TSV files that contained V, D and J usage frequencies,
and the distribution of CDR3 length,
and summarize them into HTML report (`stats_report.html`).



.. code-block:: bash
    
    pydair stats -i ./results/joinsim.vdj.pydair ./results/joinsim.sub.vdj.pydair \
                 -n whole_data subset \
                 -o ./results/stats \
                 --contain_ambiguous_D \
                 --estimate-vdj-combination










Analysis of human HIV-1-neutralizing antibodies
===============================================


We show the precedures for repertoire diversity study of
human immunoglobulin heavy (IgH) chains from B cell with PyDAIR.
The IgH sequences were sequenced from the two donors IVAI84 and N152 using 454 pyrosequencing
in `Zhu et al <http://www.pnas.org/content/110/16/6470.long>`_.
IgH sequence in IAVI84 donor is broadly contained neutralizing antibodies,
and N152 is the brodly neutralizing antibody 10E8 was recently identified in HIV-1-infected donor.


First, we use ``git clone`` command to download
the case study set that consists of
human germline genes in FASTA format from
`GitHub PyDAIR repository <https://github.com/bioinfoteam/PyDAIR>`_.


.. code-block:: bash
    
    git clone git@github.com:bioinfoteam/PyDAIR.git


The data are saved in `PyDAIR/casestudies/hiv`.
We use ``cd`` command to go to `hiv` directory.


.. code-block:: bash
    
    cd PyDAIR/casestudies/hiv


.. note:: To perform analysis from FASTQ file, one may need to install 
          `NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_
          and `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_.


Preparation
^^^^^^^^^^^

Before analysis, we create BLAST database with human
germline gene sequences using ``makeblastdb``.


.. code-block:: bash
    
    cd db
    makeblastdb -in human.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


The IgH sequencing data for the two donors are available on
`NCBI SRA <www.ncbi.nlm.nih.gov/sra>`_ with the accession number of SRR654171 and SRR654169.
We use 
`NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_
to downlaod Rep-Seq data and covert them to FASTQ format file.

.. code-block:: bash
    
    prefetch SRR654171
    prefetch SRR654169
    fastq-dump SRR654171 -O ./data/
    fastq-dump SRR654169 -O ./data/


High-throughput sequencing data generally contains low qualities reads.
We use 
`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
to removed the low quality reads.


.. code-block:: bash
    
    fastqc ./data/SRR654171.fastq -o ./data/ -q --nogroup
    fastqc ./data/SRR654169.fastq -o ./data/ -q --nogroup
    
    trimmomatic SE -phred33 ./data/SRR654171.fastq ./data/SRR654171.qc.fastq HEADCROP:10 TRAILING:20 MINLEN:100
    trimmomatic SE -phred33 ./data/SRR654169.fastq ./data/SRR654169.qc.fastq HEADCROP:10 TRAILING:20 MINLEN:100
    
    fastqc ./data/SRR654171.qc.fastq -o ./data/ -q --nogroup
    fastqc ./data/SRR654169.qc.fastq -o ./data/ -q --nogroup


After trimming of low quality bases and removing low short sequences,
we convert FASTQ format file to FASTA format file
with ``awk`` and ``sed`` commands.


.. code-block:: bash
    
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/SRR654171.fastq | sed -e 's/^@/\>/' > ./data/SRR654171.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/SRR654169.fastq | sed -e 's/^@/\>/' > ./data/SRR654169.fa





Analysis
^^^^^^^^

We use ``pydair parse`` command to assign VDJ genes and determine CDR3 sequence.


.. code-block:: bash
    
    pydair parse -s human -q ./data/SRR654169.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/SRR654171
    pydair parse -s human -q ./data/SRR654169.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/SRR654169


Then, we use ``pydair stats`` command to summarize the analysis results.
All summarized data are saved into `results` directory with prefix `stats`.


.. code-block:: bash
    
    pydair stats -i ./result/SRR654171.vdj.pydair ./result/SRR654169.vdj.pydair \
                 -n N152 TIAVI84 \
                 -o ./result/stats \
                 --contain_ambiguous_D \
                 --estimate-vdj-combination





