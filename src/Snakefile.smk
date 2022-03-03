import sys
import pandas as pd
sys.path.insert(1, "/home/joel/Documents/reusablecode")

from pop_frequeny_gnomad import ExploreGnomad
genes =[
    "PHAX", "GABRA4", "ZNF488", "ZNF773", "BRIP1", "SLC22A11", "SLC29A3",
    "WDR54", "MTHFSD", "FXYD6-FXYD2", "XRN2", "PWP2", "COL18A1", "RAD50",
    "ALOX12B", "SORD", "MAP3K15", "DAAM2", "PCNT", "SPTA1", "BAZ2B", "RANBP6",
    "ZHX3"
]

assocs = [
    "test"
]

include: "Filtering.smk"
include: "Analysis.smk"

rule all:
    input:
        expand("run_folder/results/concat/fish_burden_{assoc}.tsv", assoc = assocs)
