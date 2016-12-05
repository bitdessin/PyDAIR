============
Case studies
============


The procedures of several case studies with PyDAIR are shown in this page.
The study sets can be downloaded via :command:`git clone` command from
`GitHub PyDAIR repository <https://github.com/bioinfoteam/PyDAIR>`_.


.. code-block:: bash
    
    git clone git@github.com:bioinfoteam/PyDAIR.git


All study sets are saved in the :file:`casestudies` directory
of the downloaded :file:`PyDAIR` directory.


.. code-block:: bash
    
    cd PyDAIR
    cd casestudies
    ls
    ## hiv       mouse     simdata   zebrafish



+-------------+--------+-------------------------------------------------------------+
| dataset     | samles | description                                                 |
+=============+========+=============================================================+
| *simdata*   |      1 | Artificial sequences of human IgH.                          |
+-------------+--------+-------------------------------------------------------------+
| *hiv*       |      2 | Human HIV-1-neutralizing antibody repertoires.              |
+-------------+--------+-------------------------------------------------------------+
| *mose*      |      2 | Mouse immunoglobulin heavy chain.                           |
+-------------+--------+-------------------------------------------------------------+
| *zebrafish* |     14 | Immunoglobulin heavy chain of Zebrafish IgM and IgZ.        |
+-------------+--------+-------------------------------------------------------------+



--------------------------------------------------------------------







Analysis of `simdata` dataset
=============================


The *simdata* dataset consists of 100,000 artificial sequences of human IgH.
which were generated JoinSimulation\ [#Russ2015]_.
The raw data are saved as compressed TSV format in :file:`data/simdata.txt`.
We use :command:`bzip2` command to decompress the raw data,
and then convert the raw data into FASTA format using the custom Python script.


.. code-block:: bash
    
    cd PyDAIR/casestudies/simdata
    bzip2 -d data/simdata.txt.bz2
    
    python ./bin/convert_csv_to_fastq.py ./data/simdata.txt ./data/simdata.fa


Sequences of V, D, and J genes of human IgH are saved in the :file:`db` directory.
We use :command:`makeblast` command to create BLAST databases of V, D, and J genes.


.. code-block:: bash
    
    cd db
    makeblastdb -in human.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


After we prepared IgH sequences (:file:`simdata.fa`) and BLAST databases
(:file:`ighv` :file:`ighd` :file:`ighj`), we use :command:`pydair parse`
command to identify V, D, J, and CDR3 segments.


.. code-block:: bash
    
    pydair parse -q data/simdata.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o results/simdata


The result will be saved into :file:`results/simdata.vdj.pydair`.
In addition, the digest version of result will be saved into
:file:`results/simdata.vdj.pydair.simple` with TSV format.


Then, we use :command:`pydair stats` to summarize the anlaysis results.


.. code-block:: bash
    
    pydair stats -i ./results/simdata.vdj.pydair \
                 -n simdata \
                 -o ./results/stats \
                 --estimate-vdj-combination



The summarized results are saved into :file:`./restuls` directory with
the prefix of :file:`stats`.
The HTML report saved in :file:`./result/stats.report.html` (:download:`simdata_report.html`).



..  
    python ./bin/calc_accuracy_details.py ./data/simdata.txt \
                                  ./results/simdata.vdj.pydair \
                                  ./results/simdata.stats.p
    To evaluate the relations between the number of sequences and execution time,
    we create some subsets.
    
    head -n   2000 ./data/simdata.fa > ./data/simdata.1000.fa
    head -n  10000 ./data/simdata.fa > ./data/simdata.5000.fa
    head -n  20000 ./data/simdata.fa > ./data/simdata.10000.fa
    head -n  40000 ./data/simdata.fa > ./data/simdata.20000.fa
    head -n  80000 ./data/simdata.fa > ./data/simdata.40000.fa
    head -n 120000 ./data/simdata.fa > ./data/simdata.60000.fa
    head -n 160000 ./data/simdata.fa > ./data/simdata.80000.fa
    head -n 200000 ./data/simdata.fa > ./data/simdata.100000.fa
    
    Then, we use the custom Python script to analysis all subsetS.
    
    python ./bin/calc_exetime.py > exetime.log.txt 2>&1
    
    
    Additionally, to evaluate the relations between BLAST parameters and
    the accuracies of gene identification,
    we try six sets of BLAST parameters for V gene and six for J gene to
    analysis the first 10,0000 sequences of the original one.
    
    head -n 20000 ./data/simdata.fa  > ./data/simdata.sub.fa
    head -n 10001 ./data/simdata.txt > ./data/simdata.sub.txt
    
    python ./bin/glid_blast.py
    
    for vi in {0..5}; do
        for ji in {0..5}; do
            p=sim${vi}_${ji}
            python ./bin/calc_accuracy_details.py ./data/simdata.sub.txt \
                                          ./results/${p}.vdj.pydair \
                                          ./results/estperformance.${p}
        done
    done
    >R calc_glid_acc.R

  
--------------------------------------------------------------------






Analysis of `hiv` dataset
=========================

.. note:: To perform analysis from FASTQ file, user may need to install 
          `NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_,
          `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_,
          and `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>`_.


We show the precedures for repertoire diversity study of
human immunoglobulin heavy (IgH) chains from B cell with PyDAIR.
The IgH sequences were sequenced from the two donors IVAI84 and N152 using 454 pyrosequencing
in\ [#Zhu2013]_.
IgH sequence in IAVI84 donor is broadly contained neutralizing antibodies,
and N152 is the brodly neutralizing antibody 10E8 was recently identified in HIV-1-infected donor.


The *hiv* dataset are saved in :file:`hiv` directory.
We use :command:`cd` command to go to :file:`hiv` directory.


.. code-block:: bash
    
    cd PyDAIR/casestudies/hiv


Before analysis, we create BLAST database with human
germline gene sequences using :command:`makeblastdb`.


.. code-block:: bash
    
    cd db
    makeblastdb -in human.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


The IgH sequencing data for the two donors are available on
`NCBI SRA <www.ncbi.nlm.nih.gov/sra>`_ with the accession number of SRR654169 and SRR654171,
while SRR654169 is sequenced from IAVI84 donor
and SRR654171 is sequenced from N152 donor.
We use 
`NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_
to downlaod Rep-Seq data and covert them to FASTQ format file.


.. code-block:: bash
    
    prefetch SRR654169
    prefetch SRR654171
    fastq-dump SRR654169 -O ./data/
    fastq-dump SRR654171 -O ./data/


Both FASTQ files contain IgH and IgL sequences.
We use `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>_`
to extract the IgH sequences according to the primers.


.. code-block:: bash
    
    cutadapt -g VH15L=CCATCTCATCCCTGCGTGTCTCCGACTCAGACAGGTGCCCACTCCCAGGTGCAG \
             -g VH15L2=CCATCTCATCCCTGCGTGTCTCCGACTCAGGCAGCCACAGGTGCCCACTCC \
             -g VH124=CCATCTCATCCCTGCGTGTCTCCGACTCAGCAGCAGCTACAGGCACCCACGC \
             -g VH169=CCATCTCATCCCTGCGTGTCTCCGACTCAGGGCAGCAGCTACAGGTGTCCAGTCC \
             --discard-untrimmed -m 300 -o ./data/SRR654169.p.fastq -O 10 -e 0.3 \
             ./data/SRR654169.fastq
    
    cutadapt -g VH35L=CCATCTCATCCCTGCGTGTCTCCGACTCAGAAGGTGTCCAGTGTGARGTGCAG \
             -g VH3L1=CCATCTCATCCCTGCGTGTCTCCGACTCAGGCTATTTTAAAAGGTGTCCAATGT \
             -g VH34L1=CCATCTCATCCCTGCGTGTCTCCGACTCAGGTGGCAGCTCCCAGATGGGTCCTGTC \
             -g VH34L3=CCATCTCATCCCTGCGTGTCTCCGACTCAGGTTGCAGTTTTAAAAGGTGTCCAGTG \
             --discard-untrimmed -m 300 -o ./data/SRR654171.p.fastq -O 10 -e 0.3 \
             ./data/SRR654171.fastq
    


High-throughput sequencing data generally contains low qualities reads.
We use 
`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
to removed the low quality reads.


.. code-block:: bash
    
    fastqc ./data/SRR654169.p.fastq -o ./data/ -q --nogroup
    fastqc ./data/SRR654171.p.fastq -o ./data/ -q --nogroup
    
    trimmomatic SE -phred33 ./data/SRR654169.p.fastq ./data/SRR654169.qc.fastq TRAILING:30 MINLEN:300
    trimmomatic SE -phred33 ./data/SRR654171.p.fastq ./data/SRR654171.qc.fastq TRAILING:30 MINLEN:300
    
    fastqc ./data/SRR654169.qc.fastq -o ./data/ -q --nogroup
    fastqc ./data/SRR654171.qc.fastq -o ./data/ -q --nogroup


After trimming of low quality bases and removing low short sequences,
we convert FASTQ format file to FASTA format file
with :command:`awk` and "command:`sed` commands.


.. code-block:: bash
    
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/SRR654169.qc.fastq | sed -e 's/^@/\>/' > ./data/SRR654169.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/SRR654171.qc.fastq | sed -e 's/^@/\>/' > ./data/SRR654171.fa


As mentioned above, pydair parse was used to assign VDJ genes and define CDR3 sequences.
Analysis results were summarized via pydair stats. All the summarized data are saved into
results directory with prefix stats.
We use :command:`pydair parse` command to assign VDJ genes and determine CDR3 sequence.


.. code-block:: bash
    
    pydair parse -q ./data/SRR654169.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/SRR654169
    pydair parse -q ./data/SRR654171.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/SRR654171


Then, we use :command:`pydair stats` command to summarize the analysis results.
All summarized data are saved into :file:`results` directory with prefix `stats`.
and the summarized report were created (:download:`humanhiv_report.html`).


.. code-block:: bash
    
    pydair stats -i ./results/SRR654171.vdj.pydair ./results/SRR654169.vdj.pydair \
                 -n N152 IAVI84 \
                 -o ./results/stats \
                 --estimate-vdj-combination



--------------------------------------------------------------------



Analysis of mouse heavy chain
=============================

The datasets contains two mice of C57BL/6 and BALB/c.
Data is from\ [#Collins2015]_.

First, we use :command:`git clone` command to download
the case study set that consists of
human germline genes in FASTA format from
`GitHub PyDAIR repository <https://github.com/bioinfoteam/PyDAIR>`_.


.. code-block:: bash
    
    git clone git@github.com:bioinfoteam/PyDAIR.git


The data are saved in :file:`PyDAIR/casestudies/mouse`.
We use :command:`cd` command to go to :file:`hiv` directory.


.. code-block:: bash
    
    cd PyDAIR/casestudies/mouse


Before analysis, we create BLAST database with human
germline gene sequences using :command:`makeblastdb`.


.. code-block:: bash
    
    cd db
    makeblastdb -in mouse.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in mouse.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in mouse.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../



.. code-block:: bash
    
    wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR849/ERR849859/ERR849859.fastq.gz
    wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR849/ERR849860/ERR849860.fastq.gz
    gunzip ERR849859.fastq.gz
    gunzip ERR849860.fastq.gz


High-throughput sequencing data generally contains low qualities reads.
We use 
`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
to removed the low quality reads.


.. code-block:: bash
    
    fastqc ./data/ERR849859.fastq -o ./data/ -q --nogroup
    fastqc ./data/ERR849860.fastq -o ./data/ -q --nogroup
    
    trimmomatic SE -phred33 ./data/ERR849859.fastq ./data/ERR849859.qc.fastq HEADCROP:10 TRAILING:20 MINLEN:100
    trimmomatic SE -phred33 ./data/ERR849860.fastq ./data/ERR849860.qc.fastq HEADCROP:10 TRAILING:20 MINLEN:100
    
    fastqc ./data/ERR849859.qc.fastq -o ./data/ -q --nogroup
    fastqc ./data/ERR849860.qc.fastq -o ./data/ -q --nogroup


After trimming of low quality bases and removing low short sequences,
we convert FASTQ format file to FASTA format file
with :command:`awk` and "command:`sed` commands.


.. code-block:: bash
    
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/ERR849859.fastq | sed -e 's/^@/\>/' > ./data/ERR849859.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/ERR849860.fastq | sed -e 's/^@/\>/' > ./data/ERR849860.fa


As mentioned above, pydair parse was used to assign VDJ genes and define CDR3 sequences.
Analysis results were summarized via pydair stats. All the summarized data are saved into
results directory with prefix stats.
We use :command:`pydair parse` command to assign VDJ genes and determine CDR3 sequence.


.. code-block:: bash
    
    pydair parse -q ./data/ERR849859.fa \
                 -v ./db/mouse.ighv.fa -d ./db/mouse.ighd.fa -j ./db/mouse.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/ERR849859
    pydair parse -q ./data/ERR849860.fa \
                 -v ./db/mouse.ighv.fa -d ./db/mouse.ighd.fa -j ./db/mouse.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/ERR849860


Then, we use :command:`pydair stats` command to summarize the analysis results.
All summarized data are saved into :file:`results` directory with prefix `stats`.


.. code-block:: bash
    
    pydair stats -i ./results/ERR849859.vdj.pydair ./results/ERR849860.vdj.pydair \
                 -n ERR849859 ERR849860 \
                 -o ./results/stats \
                 --estimate-vdj-combination
    

--------------------------------------------------------------------


Analysis of zebrafish antibody repertoires
==========================================

.. note:: To perform analysis from FASTQ file, one may need to install 
          `NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_
          and `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_.

We show the precedures for repertoire diversity study of
zebrafish immunoglobulin in IgM and IgZ with PyDAIR.
The IgZ and IgM sequences were sequenced from 14 zebrafish\ [#Weinstein2009]_.

First, we used :command:`git clone` command to download the case study set that consist of
zebrafish germline genes in FASTA format from
`GitHub PyDAIR repository <https://github.com/bioinfoteam/PyDAIR>`_.


.. code-block:: bash
    
    git clone git@github.com:bioinfoteam/PyDAIR.git


The data are saved in :file:`PyDAIR/casestudies/zebrafish`.
We use :command:`cd` command to go to :file:`zebrafish` directory.


.. code-block:: bash
    
    cd PyDAIR/casestudies/zebrafish



Before analysis, we create BLAST database with human
germline gene sequences using :command:`makeblastdb`.


.. code-block:: bash
    
    cd db
    makeblastdb -in zebrafish.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in zebrafish.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in zebrafish.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


The IgH sequencing data for the two donors are available on
`NCBI SRA <www.ncbi.nlm.nih.gov/sra>`_ with the accession number of SRR654169 and SRR654171,
while SRR654169 is sequenced from IAVI84 donor
and SRR654171 is sequenced from N152 donor.
We use 
`NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_
to downlaod Rep-Seq data and covert them to FASTQ format file.

.. code-block:: bash
    
    sra=("SRR017328" "SRR017329" "SRR017330" "SRR017331" "SRR017332" "SRR017333" "SRR017334" \
         "SRR017335" "SRR017336" "SRR017337" "SRR017338" "SRR017339" "SRR017340" "SRR017341")
    
    for sid in ${sra[@]}
    do
        prefetch ${sid}
        fastq-dump ${sid} -O ./data/
    done



Both FASTQ files contain IgH and IgL sequences.
We use `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>_`
to extract the IgH sequences according to the primers.

.. code-block:: bash   
    
    for sid in ${sra[@]}
    do
        cutadapt -g IGM=TGCACTGAGACAAACCGAAG -g IGZ=TCAGAGGCCAGACATCCAAT \
                 --discard-untrimmed -m 300 -o ./data/${sid}.p.fastq -O 10 -e 0.3 \
                 --info-file ./results/${sid}.primers.info.txt \
                 ./data/${sid}.fastq


.. code-block:: bash
        
        python ./bin/read_classify.py ./results/${sid}.primers.info.txt \
                                  ./data/${sid}.fastq \
                                  ./data/${sid}.x


High-throughput sequencing data generally contains low qualities reads.
We use 
`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
to removed the low quality reads.


.. code-block:: bash
    
        trimmomatic SE -phred33 ./data/${sid}.x.igm.fq \
                   ./data/${sid}.igm.qc.fq TRAILING:30 MINLEN:100
        trimmomatic SE -phred33 ./data/${sid}.x.igz.fq \
                   ./data/${sid}.igz.qc.fq TRAILING:30 MINLEN:100
    done

.. code-block:: bash
       
    for sid in ${sra[@]}
    do
        awk 'NR % 4 == 1 || NR % 4 == 2' ./data/${sid}.igm.qc.fq | sed -e 's/^@/\>/' > ./data/${sid}.igm.fa
        awk 'NR % 4 == 1 || NR % 4 == 2' ./data/${sid}.igz.qc.fq | sed -e 's/^@/\>/' > ./data/${sid}.igz.fa
    done


After trimming of low quality bases and removing low short sequences,
we convert FASTQ format file to FASTA format file
with :command:`awk` and "command:`sed` commands.



.. code-block:: bash
    
    for sid in ${sra[@]}
    do
        pydair parse -q ./data/${sid}.igm.fa \
                 -v ./db/zebrafish.ighv.fa -d ./db/zebrafish.ighd.fa -j ./db/zebrafish.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/${sid}.igm
        pydair parse -q ./data/${sid}.igz.fa \
                 -v ./db/zebrafish.ighv.fa -d ./db/zebrafish.ighd.fa -j ./db/zebrafish.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/${sid}.igz
    done
    
    pydair stats -i ./results/SRR017328.igm.vdj.pydair ./results/SRR017329.igm.vdj.pydair \
                    ./results/SRR017330.igm.vdj.pydair ./results/SRR017331.igm.vdj.pydair \
                    ./results/SRR017332.igm.vdj.pydair ./results/SRR017333.igm.vdj.pydair \
                    ./results/SRR017334.igm.vdj.pydair ./results/SRR017335.igm.vdj.pydair \
                    ./results/SRR017336.igm.vdj.pydair ./results/SRR017337.igm.vdj.pydair \
                    ./results/SRR017338.igm.vdj.pydair ./results/SRR017339.igm.vdj.pydair \
                    ./results/SRR017340.igm.vdj.pydair ./results/SRR017341.igm.vdj.pydair \
                 -n SRR017328 SRR017329 SRR017330 SRR017331 SRR017332 SRR017333 SRR017334 \
                    SRR017335 SRR017336 SRR017337 SRR017338 SRR017339 SRR017340 SRR017341 \
                 -o ./results/stats.igm --estimate-vdj-combination
    
    pydair stats -i ./results/SRR017328.igz.vdj.pydair ./results/SRR017329.igz.vdj.pydair \
                    ./results/SRR017330.igz.vdj.pydair ./results/SRR017331.igz.vdj.pydair \
                    ./results/SRR017332.igz.vdj.pydair ./results/SRR017333.igz.vdj.pydair \
                    ./results/SRR017334.igz.vdj.pydair ./results/SRR017335.igz.vdj.pydair \
                    ./results/SRR017336.igz.vdj.pydair ./results/SRR017337.igz.vdj.pydair \
                    ./results/SRR017338.igz.vdj.pydair ./results/SRR017339.igz.vdj.pydair \
                    ./results/SRR017340.igz.vdj.pydair ./results/SRR017341.igz.vdj.pydair \
                 -n SRR017328 SRR017329 SRR017330 SRR017331 SRR017332 SRR017333 SRR017334 \
                    SRR017335 SRR017336 SRR017337 SRR017338 SRR017339 SRR017340 SRR017341 \
                 -o ./results/stats.igz --estimate-vdj-combination
    
    

The HTML reports are saved in :file:`./result/stats.igm.report.html` (:download:`zebrafish_igm_report.html`),
and :file:`./result/stats.igz.report.html` (:download:`zebrafish_igz_report.html`).





References
==========

.. [#Russ2015] Russ DE, Ho KY2, Longo NS3. HTJoinSolver: Human immunoglobulin VDJ partitioning using approximate dynamic programming constrained by conserved motifs. *BMC Bioinformatics* 2015, **16**\ :170. doi: `10.1186/s12859-015-0589-x <https://dx.doi.org/10.1186/s12859-015-0589-x>`_.
.. [#Zhu2013] Zhu J, Ofek G, Yang Y, Zhang B, Louder MK, Lu G, McKee K, Pancera M, Skinner J, Zhang Z, Parks R, Eudailey J, Lloyd KE, Blinn J, Alam SM, Haynes BF, Simek M, Burton DR, Koff WC; NISC Comparative Sequencing Program, Mullikin JC, Mascola JR, Shapiro L, Kwong PD. Mining the antibodyome for HIV-1-neutralizing antibodies with next-generation sequencing and phylogenetic pairing of heavy/light chains. *Proc Natl Acad Sci U S A.* 2013, **110**\ (16):6470-5. doi: `10.1073/pnas.1219320110 <https://dx.doi.org/10.1073/pnas.1219320110>`_.
.. [#Collins2015] Collins AM, Wang Y, Roskin KM, Marquis CP, Jackson KJ. The mouse antibody heavy chain repertoire is germline-focused and highly variable between inbred strains. *Philos Trans R Soc Lond B Biol Sci.* 2015, **370**\ (1676):20140236. dio: `10.1098/rstb.2014.0236 <https://dx.doi.org/10.1098/rstb.2014.0236>`_.
.. [#Weinstein2009] Weinstein JA, Jiang N, White RA 3rd, Fisher DS, Quake SR. High-throughput sequencing of the zebrafish antibody repertoire. *Science* 2009, **324**\ (5928):807-10. doi: `10.1126/science.1170020 <https://dx.doi.org/10.1126/science.1170020>`_.



