<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [snpit](#snpit)
  - [Installation](#installation)
    - [PyPi](#pypi)
    - [Conda](#conda)
    - [Locally](#locally)
  - [Usage](#usage)
      - [VCF input and print result to screen](#vcf-input-and-print-result-to-screen)
      - [FASTA input and write result to file](#fasta-input-and-write-result-to-file)
      - [VCF input and only use records that have PASS in the FILTER field](#vcf-input-and-only-use-records-that-have-pass-in-the-filter-field)
      - [Filtering VCF based on STATUS field](#filtering-vcf-based-on-status-field)
      - [Increase threshold for calling a lineage](#increase-threshold-for-calling-a-lineage)
    - [Full usage](#full-usage)
  - [Output format](#output-format)
  - [Contributing](#contributing)
    - [Code style](#code-style)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# snpit

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)  

Whole genome SNP based identification of members of the *Mycobacterium tuberculosis* complex. Based on code originally written by Samuel Lipworth and turned into a package by Philip Fowler and Michael Hall.

`snpit` allows rapid Mycobacterial speciation of VCF files aligned to NC000962 (H37rV) and FAST(A/Q) files.

For more information please see the article;

> Lipworth S, Jajou R, de Neeling A, et al. SNP-IT Tool for Identifying Subspecies and Associated Lineages of Mycobacterium tuberculosis Complex. Emerging Infectious Diseases. 2019;25(3):482-488. doi:[10.3201/eid2503.180894](http://dx.doi.org/10.3201/eid2503.180894).

Please email <samuel.lipworth@medsci.ox.ac.uk> with any queries.

## Installation

`snpit` requires python version 3.5 or greater.

### PyPi

```bash
# not yet setup
```

### Conda

```bash
# not yet setup
```

### Locally

There are two ways of doing this: installing to your local python packages, or in a virtual environment (recommended).

First clone the repository on your local machine and move into the directory.

```bash
git clone https://github.com/philipwfowler/snpit.git
cd snpit
```

**Virtual environment**

The instructions for installing in a virtual environment are based around using [`pipenv`](https://pipenv.readthedocs.io/en/latest/), 
which is the [python recommended](https://packaging.python.org/tutorials/managing-dependencies/#managing-dependencies) way of managing dependencies.  

```bash
# check pipenv is installed and on PATH
make init
# install snpit and dependencies
make install
# make sure it is working
make test
# activate the environment and start using
pipenv shell

```

**Without virtual environment**

*Note: We strongly encourage using a virtual environment if you are installing locally.*

```bash
python3 setup.py install --user
# make sure it is working
pytest
```

## Usage

#### VCF input and print result to screen

```bash
snpit --input in.vcf
```
*Note: You do not need to specify anything special if your file is multi-sample.*
#### FASTA input and write result to file

```bash
snpit --input in.fa --output out.tsv
```

#### VCF input and only use records that have PASS in the FILTER field

```bash
snpit -i in.vcf --filter -o out.tsv
```

#### Filtering VCF based on STATUS field

This is a custom field that has been used in some CRyPTIC pipelines. It is used as a more 
fine-grained FILTER column in that some samples may pass for a position, and others may 
not.

```bash
snpit -i in.vcf --status -o out.tsv
```

#### Increase threshold for calling a lineage

```bash
snpit -i in.vcf --threshold 95
```
The threshold is the percentage of the positions known to identify this lineage that are 
found in your sample.

### Full usage

To get the full usage/help menu for `snpit` just run 

```bash
snpit --help
```

```
usage: snpit [-h] -i INPUT [-o OUTPUT] [--threshold THRESHOLD] [--filter]
             [--status] [-v]

Whole genome SNP based identification of members of the Mycobacterium
tuberculosis complex. SNP-IT allows rapid Mycobacterial speciation of VCF
files aligned to NC000962 (H37Rv).

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the VCF or FAST(A/Q) file to read and
                        classify. File can be multi-sample and/or compressed.
  -o OUTPUT, --output OUTPUT
                        Path to output results to. Default is STDOUT (-).
  --threshold THRESHOLD
                        The percentage of snps above which a sample is
                        considered to belong to a lineage. [10.0]
  --filter              Whether to adhere to the FILTER column.
  --status              Whether to adhere to the STATUS column. This is a
                        custom field that gives more fine-grained control over
                        whether a sample passes a user-defined filtering
                        criterion, even if the record has PASS in FILTER.
  -v, --version         Show the program's version number and exit.
```

## Output format

The output file is a tab-delimited file (containing a header).

```tsv
Sample  Species Lineage Sublineage      Name    Percentage
sample1 M. tuberculosis Lineage 2       N/A     beijing 91.78
sample2 M. tuberculosis Lineage 2       N/A     beijing 97.37
sample3 M. tuberculosis Lineage 4       Haarlem haarlem 100.0
```

From left to right, the columns are:
* **Sample** - the name of the sample. This is taken from the sample column heading in the VCF or the FAST(A/Q) header.
* **Species** - Species of the call.
* **Lineage** - Lineage of the call (if Mtb.).
* **Sublineage** - Sublineage of the call (if applicable).
* **Name** - name of file in the `lib/` directory where the marker variants for this call were taken from. This also relates to the common name for the lineage in some cases.
* **Percentage** - Percentage of the call's variants found in the sample.

## Contributing

We welcome any contributions. Firstly, fork this repository and clone it locally.

Next, setup `pipenv` for the project

```bash
make init
make install
make test
```

If you wish to put in a pull request to the main repository, please write a thorough 
description of the changes you have made.

### Code style

This project uses the [black](https://github.com/psf/black) code formatter. Please ensure 
any code you wish to merge has been formatted accordingly using 

```bash
make lint
```