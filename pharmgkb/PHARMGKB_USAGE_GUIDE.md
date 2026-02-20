# PharmGKB Data Files - Usage Guide

## Overview
The `pharmgkb/` directory contains rich pharmacogenomic reference data from PharmGKB (Pharmacogene Variation Consortium). These files can significantly enhance the current analysis pipeline.

## Available Files

### 1. **variantAnnotations.zip** (Main Content)
Contains 7 TSV files with ~28,467 rows of variant annotations:

#### `var_pheno_ann.tsv` (14,312 rows)
**Purpose**: Maps variants/alleles to phenotypes and clinical outcomes
- **Key Columns**: Variant/Haplotypes, Gene, Phenotype, Significance, Direction of effect, Notes
- **Use Case**: Predict clinical phenotypes and side effects from genotypes
- **Example**: "CYP2A6 rs1801272 AA+AT → decreased risk of Tobacco Use Disorder"

#### `var_drug_ann.tsv` (12,798 rows)
**Purpose**: Maps variants/alleles to drug responses and metabolism
- **Key Columns**: Variant/Haplotypes, Gene, Drug(s), PD/PK terms, Direction of effect
- **Use Case**: Determine drug response and dosing recommendations
- **Example**: "CYP3A4*17 → decreased metabolism of nifedipine"

#### `var_fa_ann.tsv`
**Purpose**: Variant functional annotation
- **Use Case**: Validate and enhance existing allele functionality data

#### `study_parameters.tsv`
**Purpose**: Supporting study metadata for research context

### 2. **drugLabels.zip**
Contains FDA/EMA drug label annotations:

#### `drugLabels.tsv` (1,357 rows)
**Purpose**: Structured drug label information with PGx biomarkers
- **Key Columns**: PharmGKB ID, Name, Source, Biomarker Flag, Testing Level, Genes, Variants
- **Use Case**: Link analysis results to approved labeling information

#### `drugLabels.byGene.tsv`
**Purpose**: Same data indexed by gene for quick lookup

### 3. **guidelineAnnotations.json** (Multiple PA files)
**Purpose**: Clinical guideline recommendations from CPIC/DPWG
- **Use Case**: Provide expert clinical guidance based on genotypes
- **Example**: "CPIC Guideline for CYP2C19 and clopidogrel dosing"

## Integration Opportunities

### ✅ **HIGH PRIORITY ENHANCEMENTS**

#### 1. **Clinical Recommendations Engine**
```python
# Add to analysis report: actionable recommendations based on diplotype
gene_data['clinical_recommendations'] = lookup_pharmgkb_recommendations(
    gene, 
    diplotype,
    identified_drugs
)
```
**Where**: `analyze_sample_variants.py` → `_add_gene_section_to_pdf()`
**Data Source**: `var_pheno_ann.tsv`, `guidelineAnnotations.json`

#### 2. **Drug-Gene Interaction Table**
```python
# Generate table: Gene → Drugs affected → Direction of effect
affected_drugs = lookup_gene_drugs_from_pharmgkb(gene)
```
**Where**: New section in PDF report
**Data Source**: `var_drug_ann.tsv`, `drugLabels.tsv`

#### 3. **Phenotype Confidence Score**
```python
# Validate called phenotypes against PharmGKB population data
confidence = validate_phenotype_against_pharmgkb(gene, diplotype)
```
**Where**: `diplotype_phenotype.py` enrichment
**Data Source**: `var_pheno_ann.tsv`, `study_parameters.tsv`

### ✅ **MEDIUM PRIORITY ENHANCEMENTS**

#### 4. **Risk/Side Effect Prediction**
Add section to PDF showing common side effects based on genotype:
- Probability of adverse events
- Drug interactions to avoid
- Dosing adjustments needed

**Where**: PDF report under "Clinical Notes" section
**Data Source**: `var_pheno_ann.tsv` (Side effect/efficacy column)

#### 5. **Literature Support**
Link each finding to PubMed IDs:
- Show evidence level (experimental vs clinical)
- Link to PMID for researchers

**Where**: Supplementary section or JSON report
**Data Source**: All files (PMID column)

#### 6. **Drug Labeling Cross-Reference**
Add FDA/EMA label information to report:
- Which drugs have relevant PGx biomarkers
- Actionability level (Prescribing Info required? Testing required?)

**Where**: "Recommendations" section
**Data Source**: `drugLabels.tsv`

### 🔄 **IMPLEMENTATION STEPS**

#### Step 1: Create PharmGKB Data Loader
```python
# New file: scripts/pharmgkb_loader.py
class PharmGKBDataLoader:
    def load_var_pheno_annotations()
    def load_var_drug_annotations()
    def load_drug_labels()
    def lookup_recommendations(gene, diplotype)
    def get_affected_drugs(gene)
    def get_side_effects(gene, phenotype)
```

#### Step 2: Cache PharmGKB Data
- Index by gene + diplotype for O(1) lookups
- Pre-load commonly used genes during initialization
- Create bidirectional mappings: Gene↔Drug, Gene↔Phenotype, Phenotype↔SideEffect

#### Step 3: Integrate into Analysis Pipeline
- Enhance `_add_gene_section_to_pdf()` with PharmGKB recommendations
- Add new "Clinical Recommendations" section
- Add "Drug-Gene Interactions" table

#### Step 4: Enrich JSON Output
```json
{
  "gene": "CYP2C19",
  "diplotype": "*2/*3",
  "phenotype": "Poor Metabolizer",
  "pharmgkb_recommendations": [
    {
      "drug": "clopidogrel",
      "recommendation": "Avoid use; increased risk of thrombotic events",
      "evidence": "CPIC Level: A",
      "pmid": "12345678"
    }
  ],
  "affected_drugs": ["clopidogrel", "warfarin", "pantoprazole"],
  "common_side_effects": ["increased bleeding risk"]
}
```

## Quick Reference: Which File to Use When

| Task | File | Key Columns |
|------|------|------------|
| Recommend drugs to avoid | `var_drug_ann.tsv` | Direction of effect, Drug(s) |
| Predict side effects | `var_pheno_ann.tsv` | Phenotype, Side effect/efficacy/other |
| Get clinical guidelines | `guidelineAnnotations.json` | Clinical recommendations |
| Check drug label info | `drugLabels.tsv` | Genes, Biomarker Flag |
| Find population evidence | `var_pheno_ann.tsv` | Study parameters, PMID |

## Example: Adding Drug Recommendations

**Before (Current)**:
```
CYP2C19 *2/*3 → Poor Metabolizer
```

**After (With PharmGKB)**:
```
CYP2C19 *2/*3 → Poor Metabolizer
↓
CLINICAL RECOMMENDATIONS:
• AVOID: clopidogrel (increased thrombotic events - CPIC Level A)
• AVOID: pantoprazole (reduced efficacy)
• ADJUST: Warfarin (increased bleeding risk)
• SAFE: Atorvastatin, Sertraline

AFFECTED DRUGS: 15+ medications documented

SIDE EFFECTS: Increased risk of bleeding events, reduced drug efficacy
```

## Files Already Downloaded
✓ All PharmGKB zips are ready to extract and use
✓ Data is current as of 2026-01-05

## Next Steps
1. Decide which enhancements to prioritize
2. Create `pharmgkb_loader.py` to parse and index the data
3. Implement lookup functions in analysis pipeline
4. Enhance PDF/JSON report generation
5. Test with your sample data
