# Standard library imports

# Related third party imports
from Bio.SeqFeature import CompoundLocation
import pyranges as pr

from BarCodeLibrary import BarCodeLibrary
from BowtieRunner import BowtieRunner
from GenBankParser import GenBankParser
from Logger import Logger
from PySamParser import PySamParser


### Example usage ###
logger = Logger()
barcodes = BarCodeLibrary("Example_Libraries/CN-32-zmo.tsv", column="spacer")
genbank = GenBankParser("GCA_003054575.1.gb")


with BowtieRunner() as bowtie:
    bowtie.make_fasta(genbank.records)
    bowtie.make_fastq(barcodes.barcodes)
    bowtie.create_index()
    bowtie.align(num_mismatches=0, num_threads=12)

    sam = PySamParser(bowtie.sam_path)

    # TODO: This is a work in progress. I am currently exploring how to create a class.
    pam = "NGNC"
    direction = "downstream"

    targets = sam.ranges.join(genbank.ranges)

    ## compile stats ##

    # barcodes with genome targets, here the genome is called "source"
    source_targets = targets[targets.Type == "source"]
    
    logger.info(f"Found {len(source_targets)} barcodes with genome targets ...")
    
    # Create a DataFrame with unique combinations of 'Barcode', 'Chromosome', 'Start', 'End'
    unique_genomic_sites = source_targets.df.drop_duplicates(subset=['Barcode', 'Chromosome', 'Start', 'End'])

    logger.info(f"Found {len(unique_genomic_sites)} unique genomic sites ...")


    # Find barcodes that appear more than once in the unique combinations
    barcode_occurrences_in_unique_sites = unique_genomic_sites['Barcode'].value_counts()
    
    logger.info(f"Found {len(barcode_occurrences_in_unique_sites)} unique barcodes ...")
    
    barcodes_in_multiple_sites = barcode_occurrences_in_unique_sites[barcode_occurrences_in_unique_sites > 1].index

    logger.info(f"Found {len(barcodes_in_multiple_sites)} barcodes that appear in multiple unique genomic sites ...")

    multiple_site_barcodes_df = source_targets[source_targets.Barcode.isin(barcodes_in_multiple_sites)].df

    logger.info(f"Found {len(multiple_site_barcodes_df)} sites where these barcodes appear ...")    

    # Add the 'Sites' column
    multiple_site_barcodes_df['Sites'] = multiple_site_barcodes_df['Barcode'].map(barcode_occurrences_in_unique_sites)

    # Convert DataFrame back to PyRanges
    multiple_site_barcodes = pr.PyRanges(multiple_site_barcodes_df)
    
    

    def convert_pam(pam):
        return pam.replace("N", "[GACT]")

    pam_search = convert_pam(pam)

    targets.PAM = targets.df.apply(
        lambda row: genbank.get_pam_sequence(row, len(pam), direction), axis=1
    )

    targets.Offset = targets.df.apply(
        lambda row: (
            row.Start - row.Start_b if row.Strand_b == "+" else row.End_b - row.End
        ),
        axis=1,
    )

    targets.Overlap = targets.df.apply(
        lambda row: min(row.End, row.End_b) - max(row.Start, row.Start_b), axis=1
    )

    targets_with_pam = targets[targets.PAM.str.contains(pam_search, regex=True)]

    if targets_with_pam.empty:
        logger.warn(f"No targets with matching PAM: '{pam}'")

    else:
        # here we're leveraging that source contains the whole genome to get the barcode occurrences
        pam_barcode_occurrences = targets_with_pam[
            targets_with_pam.Type == "source"
        ].Barcode.value_counts()

        gene_targets_with_pam = targets_with_pam[targets_with_pam.Type == "gene"]

        logger.info(f"Taken {len(sam.ranges)} barcodes from aligner ...")
        logger.info(
            f"Found {len(targets[targets.Type == 'source'])} barcodes in the genome ..."
        )
        logger.info(
            f"Found {len(pam_barcode_occurrences)} targets with matching PAM: '{pam}'"
        )
        
        gene_targets_with_pam_in_multiple_sites = gene_targets_with_pam[
            gene_targets_with_pam.Barcode.isin(barcodes_in_multiple_sites)
        ]
        
        logger.info(f"Found {len(gene_targets_with_pam_in_multiple_sites)} gene targets with matching PAM in multiple sites ...")
        
        # return the genes with pams that are NOT in multiple sites
        gene_targets_with_pam_not_in_multiple_sites = gene_targets_with_pam[
            ~gene_targets_with_pam.Barcode.isin(barcodes_in_multiple_sites)
        ]
        
        logger.info(f"Found {len(gene_targets_with_pam_not_in_multiple_sites)} gene targets with matching PAM NOT in multiple sites ...")
        logger.info(gene_targets_with_pam_not_in_multiple_sites.df)
        
        # return sites with PAMS that do NOT have genes associated with
        sites_with_pam_no_genes = targets_with_pam[
            ~targets_with_pam.Barcode.isin(gene_targets_with_pam.Barcode)
        ]
        logger.info(f"Found {len(sites_with_pam_no_genes)} sites with PAMs that are not located in a gene ...")
        if(len(sites_with_pam_no_genes) > 0):
            logger.info(sites_with_pam_no_genes.df)
            
        # return the sites without a gene with PAMs that are NOT in multiple sites
        sites_with_pam_no_genes_not_in_multiple_sites = sites_with_pam_no_genes[
            ~sites_with_pam_no_genes.Barcode.isin(barcodes_in_multiple_sites)
        ]
        logger.info(f"Found {len(sites_with_pam_no_genes_not_in_multiple_sites)} sites with PAMs that are not located in a gene and NOT in multiple sites ...")
        if(len(sites_with_pam_no_genes_not_in_multiple_sites) > 0):
            logger.info(sites_with_pam_no_genes_not_in_multiple_sites.df)
            
        