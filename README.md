# jerboa

## Data Gathering
1. Select available rodent species from [Ensembl](https://uswest.ensembl.org/info/about/species.html) 
2. Download genome sequence from [Ensembl ftp](ftp://ftp.ensembl.org/pub/release-91/fasta/)
3. Download onthologous gene relationship from Ensembl using [bioMart](http://uswest.ensembl.org/biomart/)
4. Merge and clean onthologous gene info, get ensemble gene/protein id etc, using in house python scripts.

## Gene MSA and variants analysis:    
1. Gather protein seq and ortho’s protein seq together
2. Do Multi Seq Align (MSA) by MUSCLE
3. Call variants from MSA, using in house python scripts.
4. Using provean to evaluate variants from 3, on server
5. Filter variants based on certain provean score threshold

## Variants visualization and filtering
1. Data visualization using Dash (python package)
2. Filter variants based on PROVEAN and other info

## Useful resources:
- ENSEMBL/BioMart: genome seq, gene info, ortho gene info
- MUSCLE: do MSA
- PROVEAN: variant prediction (also combine ncbi-blast, CD-hit)
- biopython [python package]: easily deal with bio data format
- Dash [python package]: single page app for visualization
- Python data analysis packages: pandas (deal with tables), matplotlib (plot), seaborn (scientific plot)
