# PhyloMInt
Pathway-based and phylogenetically adjusted quantification of competition and cooperation between microbial species. 

Takes input assembled genomes or genome scale metabolic reconstructions (GENREs) in SBML format. Using our python implementation to calculate the metabolic complementarity index<sup>1</sup> and metabolic competition index<sup>2</sup>.

\*Currently Under Review

## Getting Started
### Prerequisites & Dependencies
Dependencies required to be installed prior to running PhyloMInt.
```
Python 3.6+
    - Pandas
    - Numpy
    - NetworkX
    - python-libsbml (https://pypi.org/project/python-libsbml/)
FragGeneScan (included)
CarveMe (requires DIAMOND aligner & IBM CPLEX optimizer)
```

### Download from Git
```
git clone https://github.com/mgtools/PhyloMint.git
```

## Usage

```
usage: PhyloMInt [-h] [-i INFILE1] [-j INFILE2] [-c CARVEME] [-d DIRR] -o
                 OUTFILE [--outdir OUTDIR] [--maxcc MAXCC]
                 [--fraggenescan FRAGGENESCAN]

PhyloMInt: Takes input genomes or genome scale metabolic models in SBML
format. Calculates metabolic complementarity and cooperation indexes. Outputs
to output file, tsv format.

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE1, --infile1 INFILE1
                        Input SBML file (XML format) organism A.
  -j INFILE2, --infile2 INFILE2
                        Input SBML file (XML format) organism B.
  -c CARVEME, --carveme CARVEME
                        (Optional) Path to fasta files. Will run FragGeneScan,
                        predict genes to be used as CarveMe output.
  -d DIRR, --dirr DIRR  (Optional) Path to directory with SBML files. Will
                        compute all pairwise comparisons & use dynamic
                        programming for computational speed-up.
  -o OUTFILE, --outfile OUTFILE
                        Output file name.
  --outdir OUTDIR       Outdir path. default = cwd
  --maxcc MAXCC         Maximum number of nodes in a strongly connected
                        component (SCC) to consider in SeedSet. (default = 5)
  --fraggenescan FRAGGENESCAN
                        path to FragGeneScan. (default:
                        lib/FragGeneScan1.30/FragGeneScan)
  --version             show program's version number and exit
```

## Example Usage

### Individual genomes 

```
PhyloMInt -i <GENRE_1.xml> -j <GENRE_2.xml> 
```

### All FASTA files in a given directory
Searches for all files with the extensions fna, fa, fq, fastq to use as input to the pipeline.
```
PhyloMInt -c <path/fasta_files> --outdir <path/outdir> -o <PhyloMInt_output.tsv>
```

### All XML files in a given directory
You can use this option if you wish to use alternative genome scale metabolic reconstruction tools or manually curated metabolic networks as input (as opposed to CarveMe).
```
PhyloMInt -d <path/SBML_files> --outdir <path/outdir> -o <PhyloMInt_output.tsv>
```

## Test Example

```
PhyloMInt -c test/genomes --outdir test/ -o test/summary.tsv
```

## Version History
v0.1.0 Initial Upload

## Sources 
1. [Levy, R., Carr, R., Kreimer, A., Freilich, S., Borenstein, E. "NetCooperate: a network-based tool for inferring host-microbe and microbe-microbe cooperation." *BMC Bioinformatics*, 2015.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0588-y)
2. [Kreimer, A., Doron-Faigenboim, A., Borenstein, E., Freilich, S. "NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species." *Bioinformatics*, 2012.](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts323)


