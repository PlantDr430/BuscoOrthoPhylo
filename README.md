# BuscoOrthoPhylo
#### Wrapper for Phylogenetic construction of BUSCO orthologs

## Summary

This is a wrapper script that takes multiple BUSCO output folders or the output folder 
from Orthofinder to concatenate, align, and contruct phylogentic tree(s) of single copy 
BUSCOs/orthologs. Please be aware of directory strucuture required for BUSCO innput (see below)

## Dependencies

1. [MAFFT](https://mafft.cbrc.jp/alignment/software/linux.html) or [MUSCLE](https://www.drive5.com/muscle/)
2. [FastTree](http://www.microbesonline.org/fasttree/#Install) or [RAxMLHPC-PTHREADS](https://github.com/stamatak/standard-RAxML)

## Usage
```
optional arguments:
  -h, --help                                show this help message and exit
  -d , --directory                          Parent directory containing directory for each species
  -i , --input                              Program used for input directory data [busco|orthofinder]
  -o , --out                                Name for output directory
  -k , --key                                Keyword to help find correct BUSCO output folder [ex. -k sordario when busco output folder is run_sordario
  -t , --type                               Choice of sequence output [prot|nuc|both] (only prot available when using orthofinder)[default: prot]
  -sub , --subsample                        Percent of genes to subsample for a preliminary analysis [ex. 0.10]
  -l , --list                               File contains list of single copy BUSCOs or Orthogroups for Phylogeny on (one per line)
  --stop                                    Stop after concatenation [ac] stop after alignment [aa]
  --resume                                  Resume after concatenation [ac] resume after alignment [aa]
  -c , --cpus                               Cores to use [default: 2]
  -a , --aligner                            Choice of sequence aligner [mafft|muscle] [default: mafft]
  -ma_a  [ ...], --mafft_args  [ ...]       Additional flags, separated by commas, within quotes, no spaces after commas, to be used in MAFFT alignment [ex. 'localpair,maxiterate 0-1000']
  -mu_i , --muscle_iter                     Maximum number of iterations for MUSCLE [default: 16]
  -mu_d, --muscle_diags                     Turn on find diagonals for MUSCLE (faster for similar sequences) [default: OFF]
  -ph , --phylogeny                         Choice of sequence output [fasttree|raxml] [default: fasttree]
  -b , --boot                               Number of bootstraps. Seriously consider reducing if using raxml. [default: 1000]
  -ft_a  [ ...], --fasttree_args  [ ...]    Additional flags, separated by commas, within quotes, no spaces after commas,to be used in fasttree run [ex. 'slow,bionj,spr 4']
  -rx_a  [ ...], --raxml_args  [ ...]       Additional flags, separated by commas, within quotes, no spaces after commas,to be used in raxml run [ex. 'k,f a,t file.nwk']
  -rx_og  [ ...], --raxml_outgroup  [ ...]  Names of outgroup species to use in raxml run [ex. 'rat,mouse,fly,human']. Names must correlate to headers in alignment file
  -rx_p_sub , --raxml_prot_sub              Protein substitution model to use [see raxmlHPC -h for all -m flag options] [default: PROTCATAUTO]
  -rx_n_sub , --raxml_nuc_sub               Nucleotide substitution model to use [see raxmlHPC -h for all -m flag options] [default: GTRCAT]
  -rx_al , --raxml_algorithm                Algorithm to use [see raxmlHPC -h for all -f flag options] [default: a]
  -rx_t, --raxml_threads                    Turn off PTHREADS [default: ON]
```


For BUSCO use please set up your -d/--directory input like so:
```
Parent_directory
    |
    |
    |
    +--Genus_species  <- This name will be used in FASTA headers and Phylogenetic trees
    |       |
    |       |
    |       +--run_xxx  <- This is your BUSCO output folder
    |
    |
    +--Genus_species
            |
            |
            +--run_xxx
```
For Orthofinder use: 

Parent_directory = Results_XXX output folder from Orthofinder

(Note: Orthofinder will only use 'prot' as Orthofinder does not output any nucleotide sequences)

(Note: Currently only compatible with Orthofinder Version 2.2 or lower. Will add 2.3 and above soon)



## Default commands
Running the script with all default settings will run these commands depending on which aligner and ML construction program you use.

MAFFT
```
mafft --threads 2 genes.fasta > alignment.fasta
```
MUSCLE
```
muscle -maxiters 16 -log logfile.log -in genes.fasta -out alignment.fasta
```
FastTree
1. Proteins
```
FastTree -seed 12345 -gamma -boot 1000 -log logfile.log alignment.fasta > tree.nwk
```
2. Nucleotides
```
FastTree -nt -seed 12345 -gamma -boot 1000 -log logfile.log alignment.fasta > tree.nwk
```
RAxML (It is highly recommended to reduce the bootstraps for RAxML with -b/--boot flag or you could be waiting a very very long time)
1. Proteins
```
raxmlHPC-PTHREADS -m PROTCATAUTO -s alignment.fasta -x 12345 -p 12345 -# 1000 -k -T 2 -f a -n prot_out
```
2. Nucleotides
```
raxmlHPC-PTHREADS -m GTRCAT -s alignment.fasta -x 12345 -p 12345 -# 1000 -k -T 2 -f a -n prot_out
```

## Output files

Output Folder

  1. species
  
      * Individual species folders containing FASTA file of genes and concatenation of genes

  2. alignment
  
      * Final FASTA concatenation of all species and alignment FASTA file(s)

  3. phylogeny
  
      * Newicks tree file(s) and other associated files
