"""
Data loading and parsing utilities for TB Dashboard.
Handles GFF3, catalogue master file, and genomic coordinates.
"""
import pandas as pd
import gffutils
import os
import sqlite3
import threading
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class GeneInfo:
    """Gene information from GFF3."""
    gene_id: str
    gene_name: str
    locus_tag: str
    start: int
    end: int
    strand: str
    product: str
    chromosome: str


@dataclass
class MutationInfo:
    """Mutation information with calculated coordinates."""
    variant: str
    gene: str
    mutation: str
    genomic_position: int
    gene_relative_position: Optional[int]  # Position relative to gene start
    gene_start: int
    gene_end: int
    strand: str
    drug: Optional[str]
    resistance_profile: Optional[str]
    tier: Optional[str]
    effect: Optional[str]
    confidence: Optional[str]


class DataLoader:
    """Load and parse TB genomic data files."""
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir
        self.gff3_path = os.path.join(data_dir, "h37rv.gff3")
        self.catalogue_path = os.path.join(data_dir, "catalogue_master_file.txt")
        self.genomic_coords_path = os.path.join(data_dir, "genomic_coordinates.txt")
        self.fasta_path = os.path.join(data_dir, "h37rv.fasta")
        self.fai_path = os.path.join(data_dir, "h37rv.fasta.fai")
        
        self.gene_db: Optional[gffutils.Database] = None
        self.genes_cache: Dict[str, GeneInfo] = {}
        self.catalogue_df: Optional[pd.DataFrame] = None
        self.genomic_coords_df: Optional[pd.DataFrame] = None
        
    def get_gene_db(self) -> gffutils.Database:
        """Get gene database, ensuring thread safety for SQLite."""
        db_path = self.gff3_path + ".db"
        
        # In Dash/Flask, multiple threads access this data.
        # SQLite connections can't be shared across threads.
        # We'll create a new connection if needed.
        if not hasattr(self, '_thread_local_db'):
            self._thread_local_db = threading.local()
        
        if not hasattr(self._thread_local_db, 'db'):
            if not os.path.exists(db_path):
                self.load_gff3()
            self._thread_local_db.db = gffutils.FeatureDB(db_path, keep_order=True)
        
        return self._thread_local_db.db

    def load_gff3(self) -> gffutils.Database:
        """Load GFF3 file and create gene database."""
        db_path = self.gff3_path + ".db"
        
        # Create database if it doesn't exist
        if not os.path.exists(db_path):
            self.gene_db = gffutils.create_db(
                self.gff3_path,
                dbfn=db_path,
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True
            )
        
        # For the current thread, initialize the DB
        return self.get_gene_db()
    
    def get_gene_info(self, gene_name: str) -> Optional[GeneInfo]:
        """Get gene information by name or locus tag."""
        if gene_name in self.genes_cache:
            return self.genes_cache[gene_name]
        
        db = self.get_gene_db()
        
        # Try searching by gene name
        try:
            genes = list(db.features_of_type('gene'))
            for gene in genes:
                name = gene.attributes.get('Name', [''])[0]
                locus = gene.attributes.get('locus_tag', [''])[0]
                
                if name == gene_name or locus == gene_name:
                    # Get product from CDS
                    product = ""
                    cds = next(db.children(gene, featuretype='CDS'), None)
                    if cds:
                        product = cds.attributes.get('product', [''])[0]
                    
                    gene_info = GeneInfo(
                        gene_id=gene.id,
                        gene_name=name,
                        locus_tag=locus,
                        start=gene.start,
                        end=gene.end,
                        strand=gene.strand,
                        product=product,
                        chromosome=gene.seqid
                    )
                    self.genes_cache[gene_name] = gene_info
                    return gene_info
        except Exception as e:
            print(f"Error finding gene {gene_name}: {e}")
        
        return None
    
    def load_catalogue(self) -> pd.DataFrame:
        """Load catalogue master file."""
        if self.catalogue_df is None:
            self.catalogue_df = pd.read_csv(
                self.catalogue_path,
                sep='\t',
                dtype=str
            )
        return self.catalogue_df
    
    def load_genomic_coordinates(self) -> pd.DataFrame:
        """Load genomic coordinates file."""
        if self.genomic_coords_df is None:
            self.genomic_coords_df = pd.read_csv(
                self.genomic_coords_path,
                sep='\t',
                dtype=str
            )
        return self.genomic_coords_df
    
    def get_drug_resistance_info(self, gene_name: str, mutation: Optional[str] = None) -> pd.DataFrame:
        """Get drug resistance information for a gene/mutation."""
        catalogue = self.load_catalogue()
        
        # Gene synonyms mapping (Catalogue name -> GFF3 names/locus_tags)
        # Verified against H37Rv GFF3
        synonyms = {
            "gyrB": ["Rv0005", "gyrB"],
            "gyrA": ["Rv0006", "gyrA"],
            "fgd1": ["Rv0407", "fgd1"],
            "rpoB": ["Rv0667", "rpoB"],
            "mmpL5": ["Rv0676c", "mmpL5"],
            "mmpS5": ["Rv0677c", "mmpS5"],
            "Rv0678": ["Rv0678"],
            "rpsL": ["Rv0682", "rpsL"],
            "fbiC": ["Rv1173", "fbiC"],
            "atpE": ["Rv1305", "atpE"],
            "rrs": ["Rvnr01", "rrs"],
            "fabG1": ["Rv1483", "fabG1"],
            "inhA": ["Rv1484", "inhA"],
            "katG": ["Rv1908c", "katG"],
            "pncA": ["Rv2043c", "pncA"],
            "eis": ["Rv2416c", "eis"],
            "ahpC": ["Rv2428", "ahpC"],
            "fbiA": ["Rv3261", "fbiA"],
            "fbiB": ["Rv3262", "fbiB"],
            "ddn": ["Rv3547", "ddn"],
            "embB": ["Rv3795", "embB"],
            "ethA": ["Rv3854c", "ethA"]
        }
        
        # Determine which catalogue 'gene' to look for
        search_genes = [gene_name.lower()]
        for cat_name, syn_list in synonyms.items():
            if gene_name.lower() in [s.lower() for s in syn_list]:
                if cat_name.lower() not in search_genes:
                    search_genes.append(cat_name.lower())
                for s in syn_list:
                    if s.lower() not in search_genes:
                        search_genes.append(s.lower())

        # Filter by gene
        mask = catalogue['gene'].str.lower().isin(search_genes)
        
        if mutation:
            # mutation might be a full variant like "Rv0678_p.Asp47fs" 
            # or just the mutation part like "p.Asp47fs"
            mutation_lower = mutation.lower()
            
            # 1. Try matching full variant in the 'variant' column
            variant_mask = mask & (catalogue['variant'].str.lower() == mutation_lower)
            
            # 2. Try matching mutation part in the 'mutation' column
            mut_part = mutation_lower
            if '_' in mutation_lower:
                mut_part = mutation_lower.split('_', 1)[1]
            
            mutation_column_mask = mask & (catalogue['mutation'].str.lower() == mut_part)
            
            # 3. Special case for LoF
            if "lof" in mutation_lower:
                lof_mask = mask & (catalogue['mutation'].str.lower() == "lof")
                final_mask = variant_mask | mutation_column_mask | lof_mask
            else:
                final_mask = variant_mask | mutation_column_mask
                
            result = catalogue[final_mask].copy()
            return result
        
        result = catalogue[mask].copy()
        return result
    
    def calculate_relative_position(self, genomic_position: int, gene_info: GeneInfo) -> int:
        """
        Calculate position relative to gene start.
        For + strand: relative_pos = genomic_pos - gene_start
        For - strand: relative_pos = gene_end - genomic_pos
        """
        if gene_info.strand == '+':
            return genomic_position - gene_info.start
        else:
            return gene_info.end - genomic_position
    
    def calculate_genomic_position(self, relative_position: int, gene_info: GeneInfo) -> int:
        """
        Calculate genomic position from gene-relative position.
        For + strand: genomic_pos = gene_start + relative_pos
        For - strand: genomic_pos = gene_end - relative_pos
        """
        if gene_info.strand == '+':
            return gene_info.start + relative_position
        else:
            return gene_info.end - relative_position
    
    def search_mutations_by_gene(self, gene_name: str) -> pd.DataFrame:
        """Search for all mutations related to a gene."""
        genomic_coords = self.load_genomic_coordinates()
        
        # Extract gene name from variant column (e.g., "dnaA_p.Asp3Ala" -> "dnaA")
        mask = genomic_coords['variant'].str.startswith(gene_name + '_')
        result = genomic_coords[mask].copy()
        
        return result
    
    def get_mutation_details(self, variant: str) -> Dict:
        """Get detailed information about a specific mutation."""
        genomic_coords = self.load_genomic_coordinates()
        catalogue = self.load_catalogue()
        
        # Find in genomic coordinates
        gc_mask = genomic_coords['variant'] == variant
        gc_result = genomic_coords[gc_mask]
        
        # Find in catalogue
        cat_mask = catalogue['mutation'] == variant
        cat_result = catalogue[cat_mask]
        
        return {
            'genomic_coordinates': gc_result.to_dict('records') if len(gc_result) > 0 else [],
            'catalogue': cat_result.to_dict('records') if len(cat_result) > 0 else []
        }
    
    def search_genes(self, query: str) -> List[GeneInfo]:
        """Search for genes by name (partial match)."""
        db = self.get_gene_db()
        
        results = []
        query_lower = query.lower()
        
        genes = list(db.features_of_type('gene'))
        for gene in genes:
            name = gene.attributes.get('Name', [''])[0]
            locus = gene.attributes.get('locus_tag', [''])[0]
            
            if query_lower in name.lower() or query_lower in locus.lower():
                product = ""
                cds = next(db.children(gene, featuretype='CDS'), None)
                if cds:
                    product = cds.attributes.get('product', [''])[0]
                
                gene_info = GeneInfo(
                    gene_id=gene.id,
                    gene_name=name,
                    locus_tag=locus,
                    start=gene.start,
                    end=gene.end,
                    strand=gene.strand,
                    product=product,
                    chromosome=gene.seqid
                )
                results.append(gene_info)
        
        return results
    
    def get_jbrowse_config(self, region: str = None, start: int = None, end: int = None) -> dict:
        """Generate JBrowse 2 configuration with relative URLs."""
        # Use relative URLs that the Flask app will serve
        fasta_url = "/data/h37rv.fasta"
        fai_url = "/data/h37rv.fasta.fai"
        gff3_url = "/data/h37rv.gff3"
        
        assembly = {
            "name": "H37Rv NC_000962.3",
            "sequence": {
                "type": "ReferenceSequenceTrack",
                "trackId": "h37rv-seq",
                "adapter": {
                    "type": "IndexedFastaAdapter",
                    "fastaLocation": {
                        "uri": fasta_url
                    },
                    "faiLocation": {
                        "uri": fai_url
                    }
                }
            }
        }
        
        tracks = [
            {
                "type": "FeatureTrack",
                "trackId": "genes",
                "name": "Genes (GFF3)",
                "assemblyNames": ["H37Rv NC_000962.3"],
                "adapter": {
                    "type": "Gff3Adapter",
                    "gffLocation": {
                        "uri": gff3_url
                    }
                }
            }
        ]
        
        return {
            "assembly": assembly,
            "tracks": tracks,
            "defaultSession": {
                "name": "H37Rv Session",
                "view": {
                    "id": "linear-genome-view",
                    "type": "LinearGenomeView",
                    "displayedRegions": [
                        {
                            "assemblyName": "H37Rv NC_000962.3",
                            "refName": region or "NC_000962.3",
                            "start": start or 1,
                            "end": end or 10000
                        }
                    ],
                    "tracks": [
                        {
                            "type": "FeatureTrack",
                            "configuration": "genes",
                            "displays": [
                                {
                                    "type": "LinearBasicDisplay",
                                    "configuration": "genes-LinearBasicDisplay"
                                }
                            ]
                        }
                    ]
                }
            }
        }
