# Drug Database Setup

## Required File

The scripts require the **full-drug-database.xml** file from DrugBank for extended drug interaction analysis.

### Setup Instructions

1. **Download the database**: 
   - Visit [DrugBank Official Database](https://www.drugbank.ca/)
   - Download the XML file (`full-drug-database.xml`)
   - Note: This file is very large (~1.8GB) and should not be committed to git

2. **Place the file**:
   ```
   /home/genwork2/enes/pharmacogenomics/full-drug-database.xml
   ```

3. **Environment Variable (Optional)**:
   You can set an environment variable to point to the database location:
   ```bash
   export DRUGBANK_XML_PATH="/path/to/full-drug-database.xml"
   ```

### File Size Note

The `full-drug-database.xml` file is excluded from version control (.gitignore) due to its large size. This file is not included in the repository to keep it manageable.

### Usage in Scripts

Scripts like `drug_interaction_integration.py` expect this file to be present at the repository root for comprehensive drug interaction analysis. If the file is not available, basic functionality will still work with reduced features.
