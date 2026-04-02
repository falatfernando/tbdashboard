# TB Dashboard - Genomic Resistance Explorer

TB Dashboard is a web-based genomic explorer for *Mycobacterium tuberculosis*, designed to bridge the gap between complex genomic data and clinical drug resistance interpretation. Built with Python Dash, it provides researchers and clinicians with an integrated platform to visualize genomic regions, calculate variant coordinates, and cross-reference mutations with the official WHO drug resistance catalogue.

## References

- **JBrowse 2**: Diesh et al, 2023. JBrowse 2: a modular genome browser with views of synteny and structural variation. *Genome Biology* 24:74. [https://doi.org/10.1186/s13059-023-02914-z](https://doi.org/10.1186/s13059-023-02914-z)
- **WHO Mutation Catalogue (2023)**: WHO catalogue of mutations in *Mycobacterium tuberculosis* complex and their association with drug resistance, second edition. [Publication](https://www.who.int/publications/i/item/9789240082410) | [GitHub Repository](https://github.com/GTB-tbsequencing/mutation-catalogue-2023/tree/main)

## Features

- **Gene Search**: Search for genes by name or locus tag (e.g., `dnaA`, `gyrA`, `rpoB`, `katG`)
- **Genomic Visualization**: View genomic regions using JBrowse integration with GFF3 annotation
- **Drug Resistance Profiles**: Display associated drugs and resistance tiers from the WHO catalogue
- **Coordinate Calculator**: Automatically convert between:
  - Genomic coordinates (absolute position on chromosome)
  - Gene-relative coordinates (c. notation, e.g., c.102G>A)
  - Amino acid positions (p. notation, e.g., p.Asp3Ala)
- **Drill-down Details**: Explore all nucleotide changes for each mutation

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup

1. **Create virtual environment** (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. **Install dependencies**:
```bash
pip install -r requirements.txt
```

3. **Verify data files**:
Ensure the following files exist in the `data/` directory:
- `h37rv.gff3` - Genome annotation file
- `h37rv.fasta` - Reference genome sequence
- `h37rv.fasta.fai` - FASTA index file
- `catalogue_master_file.txt` - WHO drug resistance catalogue
- `genomic_coordinates.txt` - Mutation coordinates

## Usage

### Running the Application

```bash
python app.py
```

The application will start on `http://localhost:8050`

### Quick Start

1. Open your browser to `http://localhost:8050`
2. Enter a gene name in the search box (e.g., `dnaA`, `gyrA`, `rpoB`)
3. Click the search button or press Enter
4. Explore the results:
   - **Gene Info**: Genomic coordinates, strand, length
   - **Genomic Visualization**: JBrowse view of the region
   - **Drug Resistance**: Associated drugs and mutation tiers
   - **Genomic Coordinates**: All nucleotide changes
   - **Coordinate Calculator**: Detailed position calculations

### Coordinate Calculation Feature

The application automatically calculates the relationship between:

- **Genomic Position**: Absolute position on the chromosome (e.g., position 8 on NC_000962.3)
- **Gene Relative Position**: Position relative to gene start (e.g., c.7 for dnaA)
- **Amino Acid Position**: Codon position in the protein (e.g., p.Asp3)

#### Example Calculation

For `dnaA_p.Asp3Ala`:
```
Gene: dnaA
  - Genomic coordinates: NC_000962.3:1-1524 (+ strand)
  - Gene start: 1
  
Mutation: p.Asp3Ala (Aspartic acid → Alanine at amino acid position 3)
  - Amino acid position: 3
  - Nucleotide position (relative): (3-1) × 3 + 1 = 7
  - Genomic position: 1 + 7 - 1 = 7
  
Result: The mutation is at genomic position 7, which is c.7 in HGVS notation
```

For genes on the **minus strand**, the calculation is reversed:
```
Gene on - strand:
  - Genomic position = Gene end - Relative position + 1
```

## Project Structure

```
tbdashboard/
├── app.py                      # Main Dash application
├── data_utils.py               # Data loading and parsing utilities
├── coordinate_calculator.py    # Coordinate conversion utilities
├── requirements.txt            # Python dependencies
├── README.md                   # This file
└── data/
    ├── h37rv.gff3             # Genome annotation
    ├── h37rv.fasta            # Reference sequence
    ├── h37rv.fasta.fai        # FASTA index
    ├── catalogue_master_file.txt  # Drug resistance catalogue
    └── genomic_coordinates.txt    # Mutation coordinates
```

## Data Sources

- **H37Rv Reference Genome**: NC_000962.3
- **GFF3 Annotation**: RefSeq annotations via NCBI
- **Drug Resistance Catalogue**: [WHO Catalogue (2023) second edition](https://www.who.int/publications/i/item/9789240082410) for *M. tuberculosis* drug resistance.

### DataLoader Class

```python
from data_utils import DataLoader

loader = DataLoader(data_dir="data")
loader.load_gff3()
loader.load_catalogue()
loader.load_genomic_coordinates()

# Get gene information
gene_info = loader.get_gene_info("dnaA")

# Search for mutations
mutations = loader.search_mutations_by_gene("dnaA")

# Get drug resistance info
resistance = loader.get_drug_resistance_info("dnaA")
```

### CoordinateCalculator Class

```python
from coordinate_calculator import CoordinateCalculator

calc = CoordinateCalculator(data_loader)

# Calculate all coordinates for a variant
result = calc.calculate_full_coordinates("dnaA_c.102G>A", "dnaA")
# Returns: {
#   'variant': 'dnaA_c.102G>A',
#   'gene': 'dnaA',
#   'gene_info': {...},
#   'mutation': 'c.102G>A',
#   'relative_nucleotide_position': 102,
#   'genomic_position': ...,
#   'amino_acid_position': ...
# }

# Get complete mutation information
mutation_info = calc.get_mutation_with_coordinates(
    variant="dnaA_c.102G>A",
    gene_name="dnaA",
    drug="Amikacin"
)
```

## Troubleshooting

### Common Issues

1. **"Gene not found" error**:
   - Check that the GFF3 file is properly loaded
   - Try using the locus tag instead of gene name (e.g., `Rv0001` for `dnaA`)

2. **JBrowse not displaying**:
   - Ensure the GFF3 file path is correct
   - Check browser console for JavaScript errors
   - Verify that dash-jbrowse is properly installed

3. **Slow performance**:
   - The catalogue file is large (~35MB), initial load may take time
   - Consider filtering data for specific drugs/genes

### Debug Mode

Run with debug output:
```bash
python app.py
# Check console for error messages
```

## Future Enhancements

- [ ] Add support for custom VCF file upload
- [ ] Implement mutation frequency visualization
- [ ] Export results to PDF/Excel
- [ ] Batch search for multiple genes/mutations

## License

This project is open source under GNU General Public License (GPL).

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Support

For questions or issues, please open an issue on the repository.
