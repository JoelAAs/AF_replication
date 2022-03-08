# Replication of hits, dominant fisher burden test

Under `data/AF_filter/variant_matrix/annotated_vars/{gene}_var_final.tsv` I have list of variants that (CHROM, POS, ALT, REF, selected INFO etc). Under `data/AF_filter/variant_matrix/gmats/{gene}_gmat.tsv` I have a variant matrix (variant per row, and samples as columns where each value in the matrix is 0-2 or missing.

The total workfolow looks like this:
![Rulegraph](/docs/figures/rules.png)

In rule `get_gnomad` in `src/filtering.smk` I append tag AF (should I use AF_popmax instead to match better?)  I add variant frequency to the Variant list (output `run_folder/annotated/{gene}_burden_gnomad.tsv`) .

The `filter_gene` filteres genes based on some critera. The `var_df.Distance <= params.distance_max` is just filtering of distance to closest exon (in gene) and something I get from the creation of the genotype matrix (to match WES better).

Previous in the creation of the multisample VCF that we use (input for creation of genotype matrices) we already filtered on callrate < 95%, therefore I didn't include the 90% callrate filtration step.

Under `src/FisherTest.R` I rewrote the test to better match my input format. Hopefully, it should match your implementation.


Cheerio, Joel


The total workflow looks like this:
![Dag](/docs/figures/dag.png)
