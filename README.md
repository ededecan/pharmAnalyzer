# pharmAnalyzer

A comprehensive pharmacogenomics analysis framework for studying drug-gene interactions, personalized medicine, and clinical decision support.

## Project Structure

- **scripts/**: Core analysis and processing scripts
- **pharmgkb/**: PharmGKB annotation data files (phenotype, drug, functional annotations)
- **source-allele-definition/**: Gene allele definition files (15 pharmacogenes)
- **source-allele-frequency/**: Population allele frequency data
- **source-allele-functionality/**: Functional activity classifications for alleles
- **source-diplotype-phenotype/**: Diplotype-to-phenotype mappings for genes
- **star_alleles.json**: Comprehensive star allele definitions
- **pharmacogenomics_updated.json**: Processed pharmacogenomics reference data

## Key Features

- Gene-to-drug interaction mapping
- Star allele and diplotype calling
- CPIC/DPWG guideline integration
- Drug-drug interaction analysis
- Risk stratification (RED/ORANGE/YELLOW/GREEN)
- Comprehensive reporting

## Setup & Configuration

### Prerequisites

- Python 3.8+
- pandas, numpy, matplotlib

### Drug Database Configuration

The scripts support integration with drug interaction databases. You need to provide your own drug database file:

**Option 1: Use Your Own Database**
- Prepare your drug database in XML or JSON format
- Update the script configuration to point to your database location
- Supported formats: DrugBank XML, custom JSON, TSV files

**Option 2: Leave Blank**
- The scripts will operate with reduced functionality for drug interactions
- All gene-based pharmacogenomics analysis will still work

**File Placement:**
If using a drug database file, place it in the repository root and update script references accordingly.

### Running Analysis Scripts

```bash
python scripts/analyze_sample_variants.py --input <vcf_file> --output <output_dir>
```

## Data Sources

- **PharmGKB**: Pharmacogenomics Knowledge Base (https://www.pharmgkb.org/)
- **CPIC**: Clinical Pharmacogenetics Implementation Consortium
- **DPWG**: Dutch Pharmacogenetics Working Group
- **PharmacoGx References**: Standard pharmacogenomics databases

## License

Please see LICENSE file for details.

## Contributing

Contributions are welcome. Please submit pull requests with detailed descriptions of changes.

## Contact

For questions or issues, please open an issue on GitHub.
