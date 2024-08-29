import os
import datetime
from scripts.resources import Resources
from snakemake.utils import min_version, validate
from snakemake.logging import logger
from snakemake.shell import shell
import pandas as pd

def import_samples():
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"]
    
    # check if sample names match file names
    not_found = []
    for sample in SAMPLES:
        r1= f"reads/{sample}_R1_001.fastq.gz"
        r2= f"reads/{sample}_R2_001.fastq.gz"
        if not os.path.isfile(r1):
            not_found.append(r1)
        if not os.path.isfile(r2):
            not_found.append(r2)
    if len(not_found) != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"ERROR: some files not found:\n{not_found}")
    
    return SAMPLES


def comparisons():
    """
    Create pairwise comparison strings from samples.csv
    """
    sample_info = pd.read_csv("config/samples.csv")
    
    if len(sample_info["genotype"].unique()) > 1 and len(sample_info["treatment"].unique()) > 1:
        # Combine genotype and treatment to get unique conditions
        sample_info["condition"] = sample_info[["genotype","treatment"]].agg('_'.join, axis=1)

        # Get reference conditions
        reference_conditions = sample_info[sample_info["reference"] == "yes"]["condition"].unique().tolist()
        
        # Get test conditions
        test_conditions = sample_info[sample_info["reference"] != "yes"]["condition"].unique().tolist()
    elif len(sample_info["genotype"].unique()) > 1 and len(sample_info["treatment"].unique()) == 1:
        # Get reference conditions
        reference_conditions = sample_info[sample_info["reference"] == "yes"]["genotype"].unique().tolist()

        # Get test conditions
        test_conditions = sample_info[sample_info["reference"] != "yes"]["genotype"].unique().tolist()
    elif len(sample_info["genotype"].unique()) == 1 and len(sample_info["treatment"].unique()) > 1:
        # Get reference conditions
        reference_conditions = sample_info[sample_info["reference"] == "yes"]["treatment"].unique().tolist()

        # Get test conditions
        test_conditions = sample_info[sample_info["reference"] != "yes"]["treatment"].unique().tolist()
    else:
        raise ValueError("Cannot create comparisons with only one treatment and one genotype...")
    
    # Create strings for comparisons
    comparisons = []
    for test in test_conditions:
        for ref in reference_conditions:
            comparisons.append(f"{test}_vs_{ref}")
    
    return comparisons


def gtf():
    if config["viral_genome"]["apply"]:
        return resources.gtf_combined
    else:
        return resources.gtf