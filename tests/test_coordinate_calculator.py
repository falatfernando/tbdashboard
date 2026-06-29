import os
import sys

# Ensure parent directory is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from data_utils import DataLoader
from coordinate_calculator import CoordinateCalculator


def test_calculate_genomic_position():
    # Initialize DataLoader with the data directory
    loader = DataLoader("data")
    
    # Initialize CoordinateCalculator with the data loader
    calculator = CoordinateCalculator(loader)
    
    # Calculate the full coordinates for variant dnaA_c.7G>A
    result = calculator.calculate_full_coordinates("dnaA_c.7G>A", "dnaA")
    
    # Assert that genomic_position equals 7
    assert result.get("genomic_position") == 7
