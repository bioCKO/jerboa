# jerboa
1. Select available rodent species from Ensembl https://uswest.ensembl.org/info/about/species.html
2. Download genome sequence from Ensembl ftp ftp://ftp.ensembl.org/pub/release-91/fasta/
3. Download onthologous gene relationship from Ensembl using bioMart http://uswest.ensembl.org/biomart/
4. Merge and clean onthologous gene info, get ensemble gene/protein id etc, using in house python scripts.

For each gene:
	5. Gather protein seq and ortho’s protein seq together
	6. Do Multi Seq Align (MSA) by MUSCLE
	7. Call variants from MSA, using in house python scripts.
	8. Using provean to evaluate variants from 7, on server
	9. Filter variants based on certain provean score threshold

Then:
10. Gather remained variants and genes do more analysis
