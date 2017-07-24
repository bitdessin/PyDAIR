============
Case studies
============


There are 3 examples to study IgH sequence diversity by using PyDAIR written in this page.
The datasets of these three examples can be downloaded from
`GitHub PyDAIR repository <https://github.com/biunit/PyDAIR>`_
via :command:`git` command.


.. code-block:: text
    
    git clone git@github.com:biunit/PyDAIR.git


Datasets of these three exmaples, named *simdata*, *hiv*, and *fugu*,
are stored in the :file:`PyDAIR/casestudies` directory.


.. code-block:: text
    
    cd PyDAIR/casestudies
    ls
    ## fugu   hiv     simdata


+-------------+---------+--------+----------------------------------------------------------------+
| dataset     | samples | reads  | description                                                    |
+=============+=========+========+================================================================+
| *simdata*   |       1 | single | A case study for evaluating PyDAIR performances with           |
|             |         |        | artificial IgH sequences.                                      |
+-------------+---------+--------+----------------------------------------------------------------+
| *hiv*       |       2 | single | Human HIV-1-neutralizing antibody repertoires.                 |
+-------------+---------+--------+----------------------------------------------------------------+
| *fugu*      |       3 | paired | Takifugu IgT and IgM sequences.                                |
+-------------+---------+--------+----------------------------------------------------------------+






Analysis of `simdata` dataset
=============================

We show the procedure to evaluate PyDAIR performances with artificial IgH sequences.
The procedure can be described as three steps:
(i) generate artificial IgH sequences with :command:`sim` mode;
(ii) analyze the artificial IgH sequences with :command:`parse` and :command:`stats` modes;
and (iii) evaluate PyDAIR performances with :command:`eval` mode.


Artificial sequence generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use :command:`sim` mode to generate 10,000 artificial IgH sequences.
:command:`sim` mode can randomly choose single of V, D, J genes from FASTA files
which are stored in the :file:`db` directory,
and merged the VDJ genes into a sequence after random mutations.


.. code-block:: text
    
    cd simdata
    pydair sim -n 10000 -o ./data/simseq.fa \
               -v ./db/human.ighv.fa \
               -d ./db/human.ighd.fa \
               -j ./db/human.ighj.fa


The artificial IgH sequences are stored in :file:`data/simseq.fa` file.
The sequence headers include the simulated conditions such as selected VDJ genes,
5' and 3' deletions and V-D and D-J insertions.


.. code-block:: bash
    
    head ./data/simseq.fa
    ## >SEQ-000000001|IGHV3-30*02|IGHD1-14*01|IGHJ2*01|V5DEL:CAGGTGCA|V3DEL:AGA|D5DEL:G|D3DEL:|J5DEL:CTACTGG|J3DEL:CTCCTCAG|VDINS:ACAC|DJINS:CT
    ## GCTGGAGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGGGGTCCATGAGTCTCTCCTGTGCAGCGTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGCCCGCCAGGCTCCAGGCAAGGGGCTGGAGAGGGTGGCATTTATACCGTATGATGGAATTATTAAATACTATGCAGACCCGTGAAGGGCCGATCACCATCTCCACAGACATTTCCAAGAACACGCTGTATCTGCAAATGAACAGCCCGAGAGCTGAGGACGCGGCTGTGTATTACTGTGCGAAACACGTATAACCGGAACCACCTTACTGCGATCTCTGGGGCCGTGGCACCCTGGTCACAGT
    ## >SEQ-000000002|IGHV3-49*04|IGHD2/OR15-2a*01|IGHJ3*01|V5DEL:GAGGTGCAGCTG|V3DEL:AGAGA|D5DEL:A|D3DEL:ATGCC|J5DEL:TGATGC|J3DEL:ACCGTCTCTTCAG|VDINS:ATAT|DJINS:TTGTTC
    ## GTGGAGTCTGGGGGAGGCTTGGTACAGCCAGGGCGGTTCCTGAGACTCTCCTGTACAGCTTCTGGATTCACCTTTGGTGATTAGCTATGAGCTGGGTCCGCCAGGCTCCAGGGAGGGGCTGGAGTGGGTAGGTTTCATGAGAAGCAAAGATTATGGTGGGACAACAGAATACGCCGCGTCTGTGAAGGCAGATTCACCATCTCAAGTGATGATTCCAAAAGCATCGCCTATTGCAAATGAACAGCCTGAAAACCGAGGAACAGCCGTGTATCACTGTACTATTTGAATATTGTAATAGTACTACTTTCTTTGTTCTTTTGATGTCTGGGGCCAAGGGACAATGGTC
    ## >SEQ-000000003|IGHV3-66*02|IGHD2-15*01|IGHJ6*01|V5DEL:GAGGTGCAGCTG|V3DEL:GAGA|D5DEL:|D3DEL:TCC|J5DEL:ATTACT|J3DEL:CCGTCTCCTCAG|VDINS:AA|DJINS:GG
    ## GTGGAGTCTGGGGAAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCGTCAGTAGCAACTACATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGTTATTTATAGAGGTGGTGGCACATACTACGCAGATCGGTGAAGGGCCGATACACCACCTCCAGAGACAATTCCAAGAACACGCTGTATCTCAAATCACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCAAAGGATATTGTAGTGGTGGTAGCTGCTACGGACTACTGCTAAGCTATGGACGTCTGGGGGCAAGGGACACGGTCA
    ## >SEQ-000000004|IGHV7-40*03|IGHD3-10*02|IGHJ4*02|V5DEL:TTTTCAAT|V3DEL:GAGAGA|D5DEL:G|D3DEL:TAAC|J5DEL:ACTA|J3DEL:CGTCTCCTCAG|VDINS:ATCCG|DJINS:GCCCTACC
    ## AGAAAAGTCATATAATCTAAGTGTCAATCCGTGGATGTTAGATAAAATATGATATATGTAAATCATGGAATACTGGCAGCCAGCATGGTATGAATTCAGTGTGTCTAGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCATCACCTACACTGGGAACCCAACATATACCAACGGCTTCACAGGACGGTTTCTATTCTCCATGGACACCTCTGTCAGCATGGCGTATCTGCAGATCAGCAGCCTAAAGGCTGAGGACACGGCCGTGTATGACTGTATATCCGTATCACTATGTTCGGGGAGTTGTTAGCCCTACCCTTTGACTACTCGGCCAGGGAACGCTGGTCAC
    ## >SEQ-000000005|IGHV3-7*03|IGHD4-4*01|IGHJ4*03|V5DEL:GAGGTGCAG|V3DEL:AGAGA|D5DEL:T|D3DEL:C|J5DEL:GCTA|J3DEL:CGTCTCCTCAG|VDINS:GAGCTGTCT|DJINS:ATTCGC
    ## CTGGTGGAGTCTGGGGGAGGCTTGGTCTAGCCTGGGGGGTCCCGAGACTCTCCTGTGCAGGCTCGGGATTCACCTTTAGTAGCTATTGGATCAGCTGGGTCCGCCAGGCTCGAGGGAAGGGGTTGGAGTGGGTGGCCAACATAAAGATAGATGGAAGTGAGAAATACTATGTGGACTCTGTGAAGGGCCGATTTACCATCTCCAGAGACAACGCCAAGAACTCACTTATCTGCAAATGAACAGCCTGAGAGCCGAGCACACGGCCGTGTATTCCTGTGCGGAGCAGTCTGACTACAGTAACTAATTCGTCTTGGACTACTGGGGCCAAGGGACCCTGGTCAC



Sequence diversity study
^^^^^^^^^^^^^^^^^^^^^^^^

We then use :command:`parse` and :command:`stats` modes to analyze the artificial IgH sequences.
Since :command:`parse` mode internally uses BLAST to identify VDJ segments,
it is required to create BLAST databases with :command:`makeblast` command.



.. code-block:: text
    
    cd db
    makeblastdb -in human.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


Then, we analyze the artificial sequences using :command:`parse` mode
and summarizing the analyzed results using :command:`stats`.


.. code-block:: text
    
    mkdir results
    
    pydair parse -q data/simseq.fa     \
                 -v ./db/human.ighv.fa \
                 -d ./db/human.ighd.fa \
                 -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb  \
                 --d-blastdb ./db/ddb  \
                 --j-blastdb ./db/jdb  \
                 -o ./results/simseq

    pydair stats -i ./results/simseq.vdj.pydair \
                 -n simdata \
                 -o ./results/simseq \
                 --estimate-vdj-combination


The summarized statistics are saved into :file:`results` directory
with the prefix of :file:`simseq`.
In addition, the summarization report is saved as HTML format file
(:download:`simseq.report.html`).

After executions of :command:`sim`, :command:`parse`, and :command:`stats` modes,
we finally use :command:`eval` mode to calculate the number of
correctly and incorrectly VDJ identifications.


.. code-block:: text
    
    pydair eval -o ./results/eval.results.txt \
                --sim-condition ./data/simseq.fa \
                --parse-result ./results/simseq.vdj.pydair
    
    cat ./results/eval.results.txt
    ##
    ## 
    ##
    ##






Analysis of `hiv` dataset
=========================


.. note:: To perform analysis from FASTQ file, user may need to install 
          `NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_,
          `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_,
          and `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>`_.


We here show an example to analyze human IgH sequences with PyDAIR.
The IgH sequences, we will use here, were sequenced from the two donors
IVAI84 and N152 using 454 pyrosequencing in Zhu et al paper\ [#Zhu2013]_.
IgH sequences in IAVI84 donor is broadly contained neutralizing antibodies,
and N152 is the brodly neutralizing antibody 10E8 was recently identified in HIV-1-infected donor.

All data can be obtained from NCBI SRA with NCBI SRA Toolkit with the
accession numbers of SRR654169 and SRR654171.


.. code-block:: text
    
    cd hiv
    
    prefetch SRR654169
    prefetch SRR654171
    fastq-dump SRR654169 -O ./data/
    fastq-dump SRR654171 -O ./data/


Since the original data consists of both IgH and IgL sequences,
we use `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>`_
to extract the IgH sequences according to the primers.


.. code-block:: text
    
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
Here we use `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
to removed such low quality reads.


.. code-block:: text
    
    trimmomatic SE -phred33 ./data/SRR654169.p.fastq ./data/SRR654169.qc.fastq TRAILING:30 MINLEN:300
    trimmomatic SE -phred33 ./data/SRR654171.p.fastq ./data/SRR654171.qc.fastq TRAILING:30 MINLEN:300
    

After trimming of low quality bases and removing short sequences,
we convert FASTQ format file to FASTA format file
with :command:`awk` and :command:`sed` commands.


.. code-block:: text
    
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/SRR654169.qc.fastq | sed -e 's/^@/\>/' > ./data/SRR654169.fa
    awk 'NR % 4 == 1 || NR % 4 == 2' ./data/SRR654171.qc.fastq | sed -e 's/^@/\>/' > ./data/SRR654171.fa


After preprocessing of Ig-Seq data,
we then prepared BLAST databases with human VDJ gene sequences.


.. code-block:: text
    
    cd db
    makeblastdb -in human.ighv.fa -out vdb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighd.fa -out ddb -dbtype nucl -parse_seqids
    makeblastdb -in human.ighj.fa -out jdb -dbtype nucl -parse_seqids
    cd ../


Finally, we use :command:`parse` mode to assign VDJ genes and determine CDR3 sequence for each FASTA file.

.. code-block:: text
    
    pydair parse -q ./data/SRR654169.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/SRR654169
    
    pydair parse -q ./data/SRR654171.fa \
                 -v ./db/human.ighv.fa -d ./db/human.ighd.fa -j ./db/human.ighj.fa \
                 --v-blastdb ./db/vdb --d-blastdb ./db/ddb --j-blastdb ./db/jdb \
                 -o ./results/SRR654171


Then, we use :command:`stats` mode to summarize the analysis results.
All summarized data are saved into :file:`results` directory with prefix :file:`hiv`,
and the summarized report (:download:`hiv.report.html`) will be saved.


.. code-block:: text
    
    pydair stats -i ./results/SRR654171.vdj.pydair ./results/SRR654169.vdj.pydair \
                 -n N152 IAVI84 \
                 -o ./results/hiv \
                 --estimate-vdj-combination








Analysis of `fugu` dataset
==========================================

.. note:: To perform analysis from FASTQ file, one may need to install 
          `NCBI SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_
          `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>`_,
          and `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_.

The protocol shows the diversity analysis of torafugu IgM and IgT sequences with PyDAIR.

First of all, we create BLAST databases of VDJ genes of torafugu.
To decrease false positives caused by BLAST, we create databases of
D and J genes for IgM and IgT, separately.
Moreover, we also create databases of V genes for each V family (i.e., V1, V2, and V3).


.. code-block:: text
    
    #git clone git@github.com:bioinfoteam/PyDAIR.git
    cd PyDAIR/casestudies/fugu

    cd db
    makeblastdb -in V1.fa -out v1db -dbtype nucl -parse_seqids
    makeblastdb -in V2.fa -out v2db -dbtype nucl -parse_seqids
    makeblastdb -in V3.fa -out v3db -dbtype nucl -parse_seqids
    makeblastdb -in Dm.fa -out dmdb -dbtype nucl -parse_seqids
    makeblastdb -in Dt.fa -out dtdb -dbtype nucl -parse_seqids
    makeblastdb -in Jm.fa -out jmdb -dbtype nucl -parse_seqids
    makeblastdb -in Jt.fa -out jtdb -dbtype nucl -parse_seqids
    cd ../


Then, we download Ig-Seq data from DDBJ SRA using :command:`wget` command.



.. code-block:: text
    
    ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/DRA004/DRA004021/......./DRA004021_1.fastq.bz2
    ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/DRA004/DRA004021/......./DRA004021_2.fastq.bz2
    ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/DRA004/DRA004021/......./DRA004022_1.fastq.bz2
    ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/DRA004/DRA004021/......./DRA004022_2.fastq.bz2
    ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/DRA004/DRA004021/......./DRA004023_1.fastq.bz2
    ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/DRA004/DRA004021/......./DRA004023_2.fastq.bz2
        
    bzip2 -d DRA004021_1.fastq.bz2
    bzip2 -d DRA004021_2.fastq.bz2
    bzip2 -d DRA004022_1.fastq.bz2
    bzip2 -d DRA004022_2.fastq.bz2
    bzip2 -d DRA004023_1.fastq.bz2
    bzip2 -d DRA004023_2.fastq.bz2
    
    mv DRA004021_1.fastq fugu1_1.fastq
    mv DRA004021_2.fastq fugu1_2.fastq
    mv DRA004022_1.fastq fugu2_1.fastq
    mv DRA004022_2.fastq fugu2_2.fastq
    mv DRA004023_1.fastq fugu3_1.fastq
    mv DRA004023_2.fastq fugu3_2.fastq
    


We sort of Ig-Seq reads based on the primers (i.e., (V1, V2, V3) x (Ct, Cm))
using `cutadapt <http://cutadapt.readthedocs.io/en/stable/index.html>`_.


.. code-block:: text

    for (( i = 1; i < 4; ++i ))
    do
        for (( j = 1; j < 3; ++j ))
        do
            cutadapt -g  nFVH1=CTGACCCAGTCTGAACCAGT   \
                     -g  nFVH2=TGAACAGTTGACACAGCCAGC  \
                     -g  nFVH3=GCCTGAAGTAAAAAGACCTGGA \
                     -g nVhCm1=CGTTCATGGTTGGAGGGTAC   \
                     -g nVhCt1=TCTGGGAAGAAGTCGAGAGC   \
                     --info-file fugu${i}_${j}.primers.info.txt \
                     --untrimmed-output /dev/null -o /dev/null  \
                     -O 10 -e 0.2 fugu${i}_${j}.fq >> log.cutadapt.txt
        done    
    done

    for (( i = 1; i < 4; ++i ))
    do
        python ../bin/read_classify.py --fq1 fugu${i}_1.fq --fq2 fugu${i}_2.fq \
                                       --log1 fugu${i}_1.primers.info.txt      \
                                       --log2 fugu${i}_2.primers.info.txt
    done

    


High-throughput sequencing data generally contains low qualities reads.
We use `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_
to removed the low quality reads.


.. code-block:: text

    cgene=("cm" "ct")
    vgene=("v1" "v2" "v3")
    
    for (( i = 1; i < 4; ++i ))
    do
        for (( c = 0; c < ${#cgene[@]}; ++c ))
        do
            for (( v = 0; v < ${#vgene[@]}; ++v ))
            do
                java -jar ../bin/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -threads 2 \
                     -trimlog log.trimmomatic.fugu${i}.${vgene[${v}]}.${cgene[${c}]}.txt      \
                     fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.fq \
                     fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.fq \
                     fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                     fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.qc_unpaired.fq \
                     fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                     fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.qc_unpaired.fq \
                     LEADING:20 TRAILING:20 MINLEN:30 >> log.qc.txt 2>&1
            done
        done
    done


Since these data are paired-end reads,
we use PEAR merge the paired-end into single reads considering the overlaps.

.. code-block:: text

    cgene=("cm" "ct")
    vgene=("v1" "v2" "v3")
    
    for (( i = 1; i < 4; ++i ))
    do
        for (( c = 0; c < ${#cgene[@]}; ++c ))
        do
            for (( v = 0; v < ${#vgene[@]}; ++v ))
            do
                ../bin/pear -f fugu${i}_1.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                            -r fugu${i}_2.${vgene[${v}]}.${cgene[${c}]}.qc.fq \
                            -o fugu${i}.${vgene[${v}]}.${cgene[${c}]}.pear.fq \
                            -j 2 -v 10 -n 300 -p 0.05 >> log.pear.txt
            done
        done
    done



Then, convert FASTQ format to FASTA format using :command:`awk` and :command:`sed` commands.


.. code-block:: text
 
    for (( i = 1; i < 4; ++i ))
    do
        awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v1.cm.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v1.cm.fa
        awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v2.cm.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v2.cm.fa
        awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v3.cm.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v3.cm.fa
        awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v1.ct.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v1.ct.fa
        awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v2.ct.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v2.ct.fa
        awk 'NR % 4 == 1 || NR % 4 == 2' fugu${i}.v3.ct.pear.fq.assembled.fastq | sed 's/^@/>/' > fugu${i}.v3.ct.fa
    done

      


After preparations, we use :command:`parse` and :command:`stats` to
analyze sequence diversities and summarize them.


.. code-block:: text
    
    cgene=("m" "t")         # cm, ct
    vgene=("1" "2" "3")     # v1, v2, v3
    
    # Identify V, D, and J genes
    for (( i = 1; i < 4; ++i ))
    do
        for (( c = 0; c < ${#cgene[@]}; ++c ))
        do
            for (( v = 0; v < ${#vgene[@]}; ++v ))
            do
                pydair parse -q fugu${i}.v${vgene[${v}]}.c${cgene[${c}]}.fa \
                             -v ../db/V${vgene[${v}]}.fa \
                             -d ../db/D${cgene[${c}]}.fa \
                             -j ../db/J${cgene[${c}]}.fa \
                             --v-blastdb ../db/v${vgene[${v}]}db --v-evalue-cutoff 1e-90 \
                             --d-blastdb ../db/d${cgene[${c}]}db \
                             --j-blastdb ../db/j${cgene[${c}]}db --j-evalue-cutoff 1e-9 \
                             -o fugu${i}.v${vgene[${v}]}.c${cgene[${c}]}
            done
        done
    done
    
    for (( i = 1; i < 4; ++i ))
    do
        for (( c = 0; c < ${#cgene[@]}; ++c ))
        do
            cat fugu${i}.v1.c${cgene[${c}]}.vdj.pydair >  fugu${i}.c${cgene[${c}]}.pydair
            cat fugu${i}.v2.c${cgene[${c}]}.vdj.pydair >> fugu${i}.c${cgene[${c}]}.pydair
            cat fugu${i}.v3.c${cgene[${c}]}.vdj.pydair >> fugu${i}.c${cgene[${c}]}.pydair
        done
    done

    for ((c = 0; c < ${#cgene[@]}; ++c))
    do
        pydair stats -i fugu1.c${cgene[${c}]}.pydair fugu2.c${cgene[${c}]}.pydair fugu3.c${cgene[${c}]}.pydair \
                     -n Fugu1 Fugu2 Fugu3 \
                     -o fugustats_c${cgene[${c}]} \
                     --estimate-vdj-combination
    done

    






References
==========

.. [#Zhu2013] Zhu J, Ofek G, Yang Y, Zhang B, Louder MK, Lu G, McKee K, Pancera M, Skinner J, Zhang Z, Parks R, Eudailey J, Lloyd KE, Blinn J, Alam SM, Haynes BF, Simek M, Burton DR, Koff WC; NISC Comparative Sequencing Program, Mullikin JC, Mascola JR, Shapiro L, Kwong PD. Mining the antibodyome for HIV-1-neutralizing antibodies with next-generation sequencing and phylogenetic pairing of heavy/light chains. *Proc Natl Acad Sci U S A.* 2013, **110**\ (16):6470-5. doi: `10.1073/pnas.1219320110 <https://dx.doi.org/10.1073/pnas.1219320110>`_.



