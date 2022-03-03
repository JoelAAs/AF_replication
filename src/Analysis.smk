
rule concat_test:
    input:
        expand("run_folder/results/{gene}_fish_burden_{{assoc}}.tsv", gene = genes)
    output:
        "run_folder/results/concat/fish_burden_{assoc}.tsv"
    shell:
        """
        head -n 1 {input[0]} >> {output}
        for tsv in {input}
        do
            tail -1 $tsv >> {output}
        done
        """


rule fisher_test:
    input:
        filtered_variants = "run_folder/filtered/{gene}_burden_gnomad.tsv",
        gmat = "data/AF_filter/variant_matrix/gmats/{gene}_gmat.tsv",
        phenofile = "data/phenotypes/phenofile_{assoc}"
    output:
        results = "run_folder/results/{gene}_fish_burden_{assoc}.tsv",
        driver = "run_folder/results/drivers/{gene}_fish_{assoc}.tsv"
    shell:
        """
        Rscript src/FisherTest.R {input.gmat} {input.filtered_variants} {input.phenofile} \
            {wildcards.gene} {output.results} {output.driver} 
        """