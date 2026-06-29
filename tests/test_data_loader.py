import os
import sys
import dataclasses

# Ensure parent directory is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from data_utils import DataLoader


def test_data_loader_gene_info():
    # Initialize DataLoader
    loader = DataLoader("data")
    
    # Get gene info for dnaA
    gene_info = loader.get_gene_info("dnaA")
    
    # Convert the dataclass to a dictionary to assert on key existence
    gene_dict = dataclasses.asdict(gene_info)
    
    # Assertions
    assert isinstance(gene_dict, dict)
    assert 'gene_name' in gene_dict
