# TB-Dashboard: An open-source web tool for visualizing WHO-concordant genomic drug resistance mutations in Mycobacterium tuberculosis for non-bioinformaticians

**Abstract**

**Background:** The rise of drug-resistant *Mycobacterium tuberculosis* demands rapid genomic interpretation. Although next-generation sequencing is increasingly accessible, translating complex genomic variants into clinical knowledge remains a bottleneck. Non-bioinformaticians, including clinical researchers and healthcare workers, require intuitive, accessible tools to interpret mutations without relying on command-line pipelines.

**Methods:** We developed TB-Dashboard, a lightweight, open-source web application designed using the Python Dash framework to democratize TB genomic analysis. The application integrates the official World Health Organization (WHO) catalogue of mutations associated with drug resistance. For visual verification of genomic regions, it embeds the JBrowse 2 genome browser via `dash-jbrowse`, enabling live track inspection. Custom coordinate translation modules resolve mappings across different genomic strands.

**Results:** TB-Dashboard provides a unified interface with five core features: (1) a fast, index-driven gene search module matching names or locus tags; (2) interactive JBrowse-powered track visualizations displaying local features; (3) comprehensive drug resistance profiling detailing confidence levels and resistance tiers from the WHO catalogue; (4) a precise coordinate calculator translating HGVS c. or p. notations to absolute genomic positions; and (5) drill-down views displaying all associated nucleotide changes.

**Conclusion:** By bridging the gap between genomic databases and clinician-friendly visualizations, TB-Dashboard facilitates local genomic surveillance. Its lightweight, dependency-light deployment profile makes it highly suitable for diagnostic laboratories and low-resource healthcare settings globally. TB-Dashboard represents a crucial step toward accessible, data-driven tuberculosis resistance monitoring.
