rule filter_gene:
    params:
        gnomad_max = 0.01,
        distance_max = 5,
    input:
        joined_variants = "run_folder/annotated/{gene}_burden_gnomad.tsv"
    output:
        filtered_variants = "run_folder/filtered/{gene}_burden_gnomad.tsv"
    run:
        var_df = pd.read_csv(input.joined_variants, sep = "\t")
        var_df = var_df[
            (var_df.Distance <= params.distance_max) &
            (var_df.AF <= params.gnomad_max) &
            (
                    (~var_df["ExonicFunc.ensGene"].isin(["synonymous_SNV", "."])) |
                    (var_df["Func.refGene"] == "splicing")
            ) &
            (var_df.VF != "0")
        ]
        var_df.to_csv(output.filtered_variants, sep="\t", index=None)



rule get_gnomad:
    params:
        gnomad = "/media/joel/Encrypted/gnomad.genomes.r2.1.1.sites.vcf.bgz",
        terms = ["AF"]
    input:
        variant_file = "data/AF_filter/variant_matrix/annotated_vars/{gene}_var_final.tsv"
    output:
        gnomad_varinats = "run_folder/annotated/{gene}_var_gnomad.tsv",
        joined_variants = "run_folder/annotated/{gene}_burden_gnomad.tsv"
    run:
        eg = ExploreGnomad(params.gnomad, input.variant_file, chrom="CHROM", id="ID", sep = ";")
        eg.match(params.terms, output.gnomad_varinats)

        burden_df = pd.read_csv(input.variant_file, sep=";")
        gnomad_df = pd.read_csv(output.gnomad_varinats, sep="\t")

        joined_df = pd.merge(
            burden_df, gnomad_df,
            on = ["CHROM", "POS", "ALT", "REF"],
            suffixes = ("_burden", "_gnomad"))
        joined_df.to_csv(output.joined_variants, sep="\t", index=None)