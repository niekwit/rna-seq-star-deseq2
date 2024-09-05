import os

class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build
    """
    # Create genome directory
    os.makedirs("resources/", exist_ok=True)
    
    def __init__(self, genome, build):
        self.genome = genome
        self.build = build
                
        # base URLs
        base_url_ens = f"https://ftp.ensembl.org/pub/release-{build}/"
                
        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"
            elif genome == "hg38":
                name = "GRCh38"
            
            # Create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
            
            # Genomes for Biomart/Enrichr
            self.ensembl = "hsapiens_gene_ensembl"
            self.enrichr = "human"
                              
        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"
            elif genome == "mm39":
                name = "GRCm39"
                
            # Create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
            
            # Genomes for Biomart/Enrichr
            self.ensembl = "mmusculus_gene_ensembl"
            self.enrichr = "mouse"
            
        elif genome == "test":
            name = "test"
            self.fasta_url = "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
            self.gtf_url = "https://github.com/niekwit/rna-seq-star-tetranscripts/raw/main/.test/Homo_sapiens.GRCh38.111_chr22.gtf.gz"
            
        self.name = name
        # Downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        
        # Path for combined file if viral genome is also selected
        self.fasta_combined = "resources/combined.fa"
        self.gtf_combined = "resources/combined.gtf"
        
    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file
        """
        return f"resources/{os.path.basename(url).replace('.gz','')}"