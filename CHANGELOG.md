# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-06-29

### Added
- **Gene Search**: Interactive lookup for *Mycobacterium tuberculosis* genes with full annotation display (locus tag, coordinates, product).
- **JBrowse Integration**: Built-in genomic browser leveraging `dash-jbrowse` to visualize gene models, genomic tracks, and custom sequence alignments.
- **Coordinate Calculator**: Utility for translating between gene-relative (HGVS c. and p. notation) and absolute genomic coordinates.
- **WHO Catalogue Integration**: Preloaded WHO tuberculosis mutation catalogue for drug resistance annotations.
- **Packaging and CI Setup**: Included standard `setup.py` packaging, unit tests, and GitHub Actions CI workflow.
