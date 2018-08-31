# snpit

Whole genome SNP based identification of members of the Mycobacterium tuberculosis complex. Based on code originally written by Samuel Lipworth and turned into a package by Philip Fowler.

SNP-IT allows rapid Mycobacterial speciation of VCF files aligned to NC000962 (H37rV).

## How to install

First clone the repository on your local machine

```   
> git clone https://github.com/philipwfowler/snpit.git
   Cloning into 'snpit'...
   remote: Counting objects: 140, done.
   remote: Compressing objects: 100% (10/10), done.
   remote: Total 140 (delta 7), reused 13 (delta 6), pack-reused 122
   Receiving objects: 100% (140/140), 2.95 MiB | 3.58 MiB/s, done.
   Resolving deltas: 100% (58/58), done.
```   

Usage python SNP-IT.py input_file output_file

Requires: csv, BioPython

The top hit against the reference SNP libraries is printed to stdout.  A complete list of absolute and % matches against all subspecies/lineages/sub-lineages is printed to the output_file.

Please email samuel.lipworth@medsci.ox.ac.uk with any queries.
