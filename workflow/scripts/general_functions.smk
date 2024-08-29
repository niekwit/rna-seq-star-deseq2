import os
import glob
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
    if paired_end:
        for sample in SAMPLES:
            r1= f"reads/{sample}_R1_001.fastq.gz"
            r2= f"reads/{sample}_R2_001.fastq.gz"
            if not os.path.isfile(r1):
                not_found.append(r1)
            if not os.path.isfile(r2):
                not_found.append(r2)
    else:
        for sample in SAMPLES:
            r1= f"reads/{sample}.fastq.gz"
            if not os.path.isfile(r1):
                not_found.append(r1)
    if len(not_found) != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"Missing files:\n{not_found}")
    
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


def paired_end():
    """
    Checks if paired-end reads are used
    """
    # Get one fastq file
    reads = glob.glob("reads/*fastq.gz")
    if len(reads) == 0:
        reads = glob.glob("reads/*fastq.gz")
    assert len(reads) != 0, "No fastq files found..."
        
    fastq = reads[0]

    # Check file extension to see if paired-end reads are used
    if fastq.endswith("_R1_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    elif fastq.endswith("_R2_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    else:
        logger.info("Single-end reads detected...")
        return False