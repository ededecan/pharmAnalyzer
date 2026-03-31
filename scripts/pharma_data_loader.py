"""
Data loader utilities for pharmacogenomics analyzer.
Centralizes CSV loading logic to reduce duplication.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, Any, Optional


class DataLoader:
    """Unified CSV data loader for pharmacogenomics files"""
    
    @staticmethod
    def safe_read_csv(csv_file: Path, dtype=None) -> Optional[pd.DataFrame]:
        """Safely read CSV file with error handling"""
        try:
            return pd.read_csv(csv_file, dtype=dtype or str, na_filter=False)
        except Exception as e:
            print(f"  Warning: Failed to parse {csv_file}: {e}")
            return None
    
    @staticmethod
    def get_safe_value(row: pd.Series, col_idx: int, default: str = "Unknown") -> str:
        """Safely extract and clean value from row at column index"""
        try:
            if col_idx >= len(row):
                return default
            value = str(row.iloc[col_idx]).strip()
            if not value or value in ['nan', 'NaN', '']:
                return default
            return value
        except (IndexError, ValueError, AttributeError):
            return default
    
    @staticmethod
    def is_empty_value(value: str) -> bool:
        """Check if value is considered empty"""
        return not value or value in ['', 'nan', 'NaN', '-']

    @staticmethod
    def load_coordinate_mapping(tsv_file: Path) -> Dict[str, Dict[str, Dict[str, str]]]:
        """Loads master_allele_mapping.tsv and extracts exact rsIDs and Variant Names"""
        import csv
        allele_dict = {}
        
        if not tsv_file.exists():
            print(f"  Warning: Master allele mapping TSV not found at {tsv_file}")
            return allele_dict
            
        with open(tsv_file, mode='r', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                gene = row.get("Gene", "").strip()
                if not gene: continue
                
                chrom = row["Chromosome"]
                pos = row["Position_GRCh38"]
                ref_alt = row["Genomic_Ref_Alt"]
                rsid = row.get("rsID", "").strip()
                variant_name = row.get("Variant_Name", "").strip()
                
                # Extract the ALT allele (e.g., from "G > A", we want "A")
                alt_allele = ref_alt.split('>')[-1].strip() if '>' in ref_alt else ref_alt.strip()
                lookup_key = f"{chrom}_{pos}_{alt_allele}"
                
                if gene not in allele_dict:
                    allele_dict[gene] = {}
                    
                allele_dict[gene][lookup_key] = {
                    'rsid': rsid,
                    'variant_name': variant_name
                }
                
        return allele_dict


class AlleleDefinitionLoader:
    """Specialized loader for allele definition CSVs"""
    
    @staticmethod
    def load(def_dir: Path) -> Dict[str, Dict[str, Any]]:
        """Load and normalize allele definitions from CSV files"""
        allele_defs = {}
        
        for csv_file in sorted(def_dir.glob("*_allele_definition.csv")):
            gene = csv_file.stem.replace("_allele_definition", "")
            df = DataLoader.safe_read_csv(csv_file)
            
            if df is None or len(df) <= 5:
                continue
            
            try:
                # Extract rsIDs from row 3
                rsid_row = df.iloc[3]
                rsid_to_col = AlleleDefinitionLoader._extract_rsid_mapping(rsid_row)
                
                # Extract alleles starting from row 5
                alleles = AlleleDefinitionLoader._extract_alleles(df, rsid_to_col)
                
                allele_defs[gene] = {
                    'rsids': list(rsid_to_col.keys()),
                    'alleles': alleles
                }
            except Exception:
                continue
        
        return allele_defs
    
    @staticmethod
    def _extract_rsid_mapping(rsid_row: pd.Series) -> Dict[str, int]:
        """Extract rsID to column mapping from rsID row"""
        mapping = {}
        for col_idx in range(1, len(rsid_row)):
            rsid_val = DataLoader.get_safe_value(rsid_row, col_idx, default="")
            if rsid_val and not DataLoader.is_empty_value(rsid_val):
                mapping[rsid_val] = col_idx
        return mapping
    
    @staticmethod
    def _extract_alleles(df: pd.DataFrame, rsid_to_col: Dict[str, int]) -> Dict[str, Dict[str, str]]:
        """Extract allele data starting from row 5"""
        alleles = {}
        
        for allele_idx in range(5, len(df)):
            allele_name = DataLoader.get_safe_value(df.iloc[allele_idx], 0, default="")
            
            if not allele_name:
                continue
            
            allele_values = {}
            for rsid, col_idx in rsid_to_col.items():
                value = DataLoader.get_safe_value(df.iloc[allele_idx], col_idx, default="")
                if value and not DataLoader.is_empty_value(value):
                    allele_values[rsid] = value
            
            alleles[allele_name] = allele_values
        
        return alleles


class AlleleFunctionalityLoader:
    """Specialized loader for allele functionality CSVs"""
    
    @staticmethod
    def load(func_dir: Path) -> Dict[str, Dict[str, Dict[str, str]]]:
        """Load allele functionality data"""
        allele_funcs = {}
        
        for csv_file in sorted(func_dir.glob("*_allele_functionality.csv")):
            gene = csv_file.stem.replace("_allele_functionality", "")
            df = DataLoader.safe_read_csv(csv_file)
            
            if df is None:
                continue
            
            allele_funcs[gene] = AlleleFunctionalityLoader._parse_functionality(df)
        
        return allele_funcs
    
    @staticmethod
    def _parse_functionality(df: pd.DataFrame) -> Dict[str, Dict[str, str]]:
        """Parse functionality data from dataframe"""
        allele_funcs = {}
        
        for idx, row in df.iterrows():
            allele_name = DataLoader.get_safe_value(row, 0, default="")
            
            if not allele_name or not allele_name.startswith('*'):
                continue
            
            allele_funcs[allele_name] = {
                'activity': DataLoader.get_safe_value(row, 1),
                'biochemical_function': DataLoader.get_safe_value(row, 2),
                'clinical_function': DataLoader.get_safe_value(row, 3),
                'summary': DataLoader.get_safe_value(row, 7, default=""),
                'strength_of_evidence': DataLoader.get_safe_value(row, 6, default="")
            }
        
        return allele_funcs


class DiplotypePhenotypeLoader:
    """Specialized loader for diplotype-phenotype CSVs"""
    
    @staticmethod
    def load(dipl_dir: Path) -> Dict[str, Dict[str, Dict[str, str]]]:
        """Load diplotype-phenotype mapping"""
        diplotype_phenos = {}
        
        for csv_file in sorted(dipl_dir.glob("*_diplotype_phenotype.csv")):
            gene = csv_file.stem.replace("_diplotype_phenotype", "")
            df = DataLoader.safe_read_csv(csv_file)
            
            if df is None:
                continue
            
            diplotype_phenos[gene] = DiplotypePhenotypeLoader._parse_phenotypes(df)
        
        return diplotype_phenos
    
    @staticmethod
    def _parse_phenotypes(df: pd.DataFrame) -> Dict[str, Dict[str, str]]:
        """Parse phenotype data from dataframe"""
        phenotypes = {}
        
        for idx, row in df.iterrows():
            try:
                diplotype = DataLoader.get_safe_value(row, 0, default="")
                
                if not diplotype or '/' not in diplotype:
                    continue
                
                phenotypes[diplotype] = {
                    'phenotype': DataLoader.get_safe_value(row, 2),
                    'ehr_priority': DataLoader.get_safe_value(row, 3)
                }
            except (IndexError, ValueError):
                continue
        
        return phenotypes


class AlleleFrequencyLoader:
    """Specialized loader for allele frequency CSVs"""
    
    @staticmethod
    def load(freq_dir: Path) -> Dict[str, Dict[str, str]]:
        """Load allele frequencies"""
        allele_freqs = {}
        
        for csv_file in sorted(freq_dir.glob("*_allele_frequency.csv")):
            gene = csv_file.stem.replace("_allele_frequency", "")
            df = DataLoader.safe_read_csv(csv_file)
            
            if df is None:
                continue
            
            allele_freqs[gene] = AlleleFrequencyLoader._parse_frequencies(df)
        
        return allele_freqs
    
    @staticmethod
    def _parse_frequencies(df: pd.DataFrame) -> Dict[str, str]:
        """Parse frequency data from dataframe"""
        frequencies = {}
        
        for idx, row in df.iterrows():
            try:
                allele_name = DataLoader.get_safe_value(row, 0, default="")
                
                if not allele_name:
                    continue
                
                frequency = DataLoader.get_safe_value(row, 1, default="N/A")
                frequencies[allele_name] = frequency
            except (IndexError, ValueError):
                continue
        
        return frequencies