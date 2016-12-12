# CSI

Complementary Sequence Index (CSI) is pipeline to quantitate RNA pairing
capacity of orientation-opposite complementary sequences across circRNA-flanking
introns.

## Schema
![pipeline](https://github.com/YangLab/CSI/blob/master/schema.png)

## Requirements
* [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279671/#introduction.Source_tarball) (>=2.2.30+)
* [pysam](http://pysam.readthedocs.org/en/latest/) (>=0.8.4pre)

## Installation
```bash
git clone https://github.com/YangLab/CSI
python ./setup.py install
```

## Usage
```
usage: CSI [-h] [-v] -g GENOME [-l LENGTH] [-p THREAD] [-o OUTPUT] [--tmp]
           circ_file

positional arguments:
  circ_file             The circular RNA file (CIRCexplorer format)

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -g GENOME, --genome GENOME
                        Genome fasta file.
  -l LENGTH, --length LENGTH
                        Minimum pair length. [default: 50]
  -p THREAD, --thread THREAD
                        Running threads. [default: 10]
  -o OUTPUT, --ouput OUTPUT
                        Output file. [default: circ_cs]
  --tmp                 Keep temporary BLAST results.
```

### Example
```bash
CSI -g hg19.fa test_circ.txt -o test_circ_csi -p 10
```

### Note
* circ_file is in the format of output of [CIRcexplorer](https://raw.githubusercontent.com/YangLab/CIRCexplorer)/[CIRCexplorer2](https://raw.githubusercontent.com/YangLab/CIRCexplorer2).
* hg19.fa is genome sequence in FASTA format.

### Output
* test_circ_csi.txt

| Field       | Description                           |
| :---------: | :------------------------------------ |
| chrom       | Chromosome                            |
| start       | Start of junction                     |
| end         | End of junction                       |
| CSI         | The max CSI score among CSs of the circRNA|
| left region | The left pair region                  |
| right region| The right pair region                 |

## Citation
**Dong R\*, Ma XK\*, Chen LL# and Yang L#. Increased complexity of circRNA expression during species evolution. RNA Biol, 2016 (Accepted)**


## License
Copyright (C) 2016 YangLab.
See the [LICENSE](https://github.com/YangLab/CSI/blob/master/LICENSE)
file for license rights and limitations (MIT).
