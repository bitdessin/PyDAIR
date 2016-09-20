[![Build Status](https://travis-ci.org/jqsunac/PyDAIR.svg?branch=master)](https://travis-ci.org/jqsunac/PyDAIR)


# PyDAIR

PyDAIR is a Python package that aims to study immunoglobulin heavy (IgH) chain diversity
based on repertoire-sequencing (Rep-Seq) data.
PyDAIR identifies the germline variable (V), diversity (D), and joining (J) genes that
used by each IgH sequence.
BLAST is used for aligning sequences to a database of known germline VDJ genes to assign VDJ.
PyDAIR supports all features as long as the two motifs, YYC and WGxG,
that located at the end of V gene and the start of J gene are know.
PyDAIR is available under the terms of the GNU license.





## INSTALLATION

PyDAIR requires Python 2.7 or Python 3.4 together with [NumPy](http://www.numpy.org/),
[Pandas](http://pandas.pydata.org/), [matplotlib](http://matplotlib.org/),
and [BioPython](http://biopython.org/) packages.
Further, PyDAIR requires [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
for aligning IgH sequence to germline databases.
PyDAIR is avaliable on the [PyPI](https://pypi.python.org/pypi/PyDAIR) repository,
as well as can be installed like any other Python package using `pip` command.

``` bash
pip install numpy --user
pip install pandas --user
pip install matplotlib --user
pip install biopython --user
pip install pydair --user
```

Installtion instructions for [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
are available on [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK279671/).
User should follow the instruction to install [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/).





## Usage

PyDAIR has two main commands that are `pydair-parseseq` and `pydair-analysis`.

| Command         | Function |
|-----------------|----------|
| pydair-parseseq | Identification of V, D and J genes that used by each IgH sequence. |
| pydair-analysis | Aggregation of the frequencies of usage of V, D and J genes, as well as extraction of CDR-H3 sequences. |


`pydair-parseseq` identifies V, D, and J genes from IgH each sequence by aligning IgH sequence to
germline (V, D, and J) database using NCBI BLAST+.
It requires IgH sequences, germline sequences, BLAST databases of germiline sequences,
and BLAST parameters. The sequences should be given by FASTA format.


```bash
 pydair-parseseq -q input_igh_sequences.fa \
                 -v v.fa                   \
                 -d d.fa                   \
                 -j j.fa                   \
                 --v-blastdb blastdb_v     \
                 --d-blastdb blastdb_d     \
                 --j-blastdb blastdb_j     \
                 -o output1
```

PyDAIR generates several files to save the intermediate results,
such as BLAST results, region that cannot be aligned to V and J genes.
The final result is saved into `output1.pydair` file.
If there several samples, `pydair-parseseq` should be run several times for each sample.

The statistical summaries are calculated by `pydair-analysis` command.

```bash
pydair-analysis -i output1.pydair output2.pydair output3.pydair  \
                -n Fugu1 Fugu2 Fugu3                             \
                -o stats_result                                  \
                --contain_ambiguous_D
```





