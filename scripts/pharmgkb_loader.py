#!/usr/bin/env python3
"""
PharmGKB Data Loader
Parses and indexes PharmGKB annotation files for fast lookup
Handles sparse/null field gracefully and extracts clinically relevant data
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Set
from collections import defaultdict


class PharmGKBDataLoader:
    """Load and index PharmGKB annotation data"""
    
    def __init__(self, pharmgkb_dir: str):
        """Initialize loader and load all PharmGKB data
        
        Args:
            pharmgkb_dir: Path to pharmgkb/ directory containing extracted TSV files
        """
        self.pharmgkb_dir = Path(pharmgkb_dir)
        
        # Initialize lookup indices
        self.gene_allele_to_phenotypes = defaultdict(list)  # (gene, allele) -> [phenotype info]
        self.gene_allele_to_drugs = defaultdict(list)       # (gene, allele) -> [drug info]
        self.gene_allele_to_functional = defaultdict(list)  # (gene, allele) -> [functional info]
        self.gene_to_drugs = defaultdict(set)               # gene -> set of drugs
        self.gene_to_side_effects = defaultdict(list)       # gene -> [side effect info]
        self.rsid_to_data = defaultdict(lambda: {'phenotypes': [], 'drugs': [], 'functional': []})  # rsID -> data
        
        self._load_phenotype_annotations()
        self._load_drug_annotations()
        self._load_functional_annotations()
        
        print(f"✓ Loaded PharmGKB data")
        print(f"  - {len(self.gene_allele_to_phenotypes)} allele-phenotype mappings")
        print(f"  - {len(self.gene_allele_to_drugs)} allele-drug mappings")
        print(f"  - {len(self.gene_allele_to_functional)} allele-functional mappings")
        print(f"  - {len(self.gene_to_drugs)} genes with drug interactions")
    
    def _load_phenotype_annotations(self) -> None:
        """Load and index phenotype annotations"""
        pheno_file = self.pharmgkb_dir / "var_pheno_ann.tsv"
        
        if not pheno_file.exists():
            print(f"⚠ Warning: Phenotype annotations not found at {pheno_file}")
            return
        
        try:
            # Load with error handling for sparse data
            df = pd.read_csv(pheno_file, sep='\t', dtype=str, na_filter=False)
            
            for _, row in df.iterrows():
                gene = str(row.get('Gene', '')).strip()
                # Extract allele/variant identifiers
                allele_raw = str(row.get('Alleles', '')).strip() or str(row.get('Variant/Haplotypes', '')).strip()
                rsid_raw = str(row.get('Variant/Haplotypes', '')).strip()
                
                # Filter allele to only include actual rsIDs (remove COSV, CM, etc.)
                # Split by '&' or ';' and keep only rs* entries
                allele_parts = [p.strip() for p in allele_raw.replace('&', ';').split(';')]
                allele = ';'.join([p for p in allele_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in allele_parts) else allele_raw
                
                # Filter rsID similarly
                rsid_parts = [p.strip() for p in rsid_raw.replace('&', ';').split(';')]
                rsid = ';'.join([p for p in rsid_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in rsid_parts) else ''
                
                drug = str(row.get('Drug(s)', '')).strip()
                phenotype = str(row.get('Phenotype', '')).strip()
                side_effect = str(row.get('Side effect/efficacy/other', '')).strip()
                direction = str(row.get('Direction of effect', '')).strip()
                significance = str(row.get('Significance', '')).strip()
                pmid = str(row.get('PMID', '')).strip()
                notes = str(row.get('Notes', '')).strip()
                sentence = str(row.get('Sentence', '')).strip()
                category = str(row.get('Phenotype Category', '')).strip()
                
                # Only keep significant or unstated findings (skip explicit 'no')
                if significance.lower() == 'no':
                    continue  # Skip only explicit negative findings
                
                # Extract side effects and link to phenotype
                if side_effect and side_effect not in ['', 'nan', 'NaN']:
                    phenotype_info = {
                        'allele': allele,
                        'gene': gene,
                        'phenotype': phenotype,
                        'side_effect': side_effect,
                        'direction': direction,
                        'drug': drug,
                        'pmid': pmid,
                        'notes': notes,
                        'sentence': sentence,
                        'category': category
                    }
                    
                    if gene and allele:
                        self.gene_allele_to_phenotypes[(gene, allele)].append(phenotype_info)
                        self.gene_to_side_effects[gene].append({
                            'allele': allele,
                            'side_effect': side_effect,
                            'direction': direction,
                            'phenotype': phenotype
                        })
                        # Also index by rsID if available
                        if rsid and rsid.startswith('rs'):
                            self.rsid_to_data[(gene, rsid)]['phenotypes'].append(phenotype_info)
        
        except Exception as e:
            print(f"⚠ Error loading phenotype annotations: {e}")
    
    def _load_drug_annotations(self) -> None:
        """Load and index drug annotations"""
        drug_file = self.pharmgkb_dir / "var_drug_ann.tsv"
        
        if not drug_file.exists():
            print(f"⚠ Warning: Drug annotations not found at {drug_file}")
            return
        
        try:
            df = pd.read_csv(drug_file, sep='\t', dtype=str, na_filter=False)
            
            for _, row in df.iterrows():
                gene = str(row.get('Gene', '')).strip()
                
                # Extract and filter allele/variant identifiers
                allele_raw = str(row.get('Alleles', '')).strip() or str(row.get('Variant/Haplotypes', '')).strip()
                rsid_raw = str(row.get('Variant/Haplotypes', '')).strip()
                
                # Filter to only include actual rsIDs (remove COSV, CM, etc.)
                allele_parts = [p.strip() for p in allele_raw.replace('&', ';').split(';')]
                allele = ';'.join([p for p in allele_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in allele_parts) else allele_raw
                
                rsid_parts = [p.strip() for p in rsid_raw.replace('&', ';').split(';')]
                rsid = ';'.join([p for p in rsid_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in rsid_parts) else ''
                
                drug = str(row.get('Drug(s)', '')).strip()
                pk_pd = str(row.get('PD/PK terms', '')).strip()
                direction = str(row.get('Direction of effect', '')).strip()
                significance = str(row.get('Significance', '')).strip()
                pmid = str(row.get('PMID', '')).strip()
                notes = str(row.get('Notes', '')).strip()
                sentence = str(row.get('Sentence', '')).strip()
                phenotype_category = str(row.get('Phenotype Category', '')).strip()
                
                # Only keep significant or unstated findings (skip explicit 'no')
                if significance.lower() == 'no':
                    continue  # Skip only explicit negative findings
                
                if gene and drug and allele:
                    # Split comma/semicolon-separated drugs to avoid duplicates
                    drug_names = [d.strip() for d in drug.replace(';', ',').split(',') if d.strip()]
                    
                    # Create separate drug_info entries for each drug
                    for drug_name in drug_names:
                        drug_info = {
                            'allele': allele,
                            'gene': gene,
                            'drug': drug_name,  # Use individual drug name
                            'direction': direction,
                            'pk_pd': pk_pd,
                            'pmid': pmid,
                            'notes': notes,
                            'sentence': sentence,
                            'category': phenotype_category
                        }
                        
                        self.gene_allele_to_drugs[(gene, allele)].append(drug_info)
                        self.gene_to_drugs[gene].add(drug_name)
                        
                        # Also index by rsID if available
                        if rsid and rsid.startswith('rs'):
                            self.rsid_to_data[(gene, rsid)]['drugs'].append(drug_info)
        
        except Exception as e:
            print(f"⚠ Error loading drug annotations: {e}")
    
    def _load_functional_annotations(self) -> None:
        """Load and index functional annotations"""
        func_file = self.pharmgkb_dir / "var_fa_ann.tsv"
        
        if not func_file.exists():
            print(f"⚠ Warning: Functional annotations not found at {func_file}")
            return
        
        try:
            df = pd.read_csv(func_file, sep='\t', dtype=str, na_filter=False)
            
            for _, row in df.iterrows():
                gene = str(row.get('Gene', '')).strip()
                
                # Extract and filter allele/variant identifiers
                allele_raw = str(row.get('Alleles', '')).strip() or str(row.get('Variant/Haplotypes', '')).strip()
                rsid_raw = str(row.get('Variant/Haplotypes', '')).strip()
                
                # Filter to only include actual rsIDs (remove COSV, CM, etc.)
                allele_parts = [p.strip() for p in allele_raw.replace('&', ';').split(';')]
                allele = ';'.join([p for p in allele_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in allele_parts) else allele_raw
                
                rsid_parts = [p.strip() for p in rsid_raw.replace('&', ';').split(';')]
                rsid = ';'.join([p for p in rsid_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in rsid_parts) else ''
                
                drug = str(row.get('Drug(s)', '')).strip()
                functional_terms = str(row.get('Functional terms', '')).strip()
                direction = str(row.get('Direction of effect', '')).strip()
                notes = str(row.get('Notes', '')).strip()
                sentence = str(row.get('Sentence', '')).strip()
                category = str(row.get('Phenotype Category', '')).strip()
                pmid = str(row.get('PMID', '')).strip()
                
                if gene and allele and functional_terms:
                    func_info = {
                        'allele': allele,
                        'gene': gene,
                        'drug': drug,
                        'functional_terms': functional_terms,
                        'direction': direction,
                        'pmid': pmid,
                        'notes': notes,
                        'sentence': sentence,
                        'category': category,
                        'type': 'functional'
                    }
                    
                    self.gene_allele_to_functional[(gene, allele)].append(func_info)
                    # Also index by rsID if available
                    if rsid and rsid.startswith('rs'):
                        self.rsid_to_data[(gene, rsid)]['functional'].append(func_info)
        
        except Exception as e:
            print(f"⚠ Error loading functional annotations: {e}")
    
    def diplotype_exists_in_db(self, gene: str, diplotype: str) -> bool:
        """Check if exact diplotype exists in PharmGKB database
        
        Args:
            gene: Gene name
            diplotype: Diplotype string (e.g., '*2/*18')
        
        Returns:
            True if exact diplotype found in database
        """
        # Check all data sources for this exact diplotype
        for key in self.gene_allele_to_drugs.keys():
            if key[0] == gene and diplotype in key[1]:
                return True
        for key in self.gene_allele_to_phenotypes.keys():
            if key[0] == gene and diplotype in key[1]:
                return True
        return False
    
    def get_drug_recommendations(self, gene: str, alleles: List[str]) -> List[Dict[str, Any]]:
        """Get drug recommendations for a gene and alleles
        
        Args:
            gene: Gene name (e.g., 'CYP2C19')
            alleles: List of alleles (e.g., ['*2', '*3'])
        
        Returns:
            List of drug recommendations with direction and evidence
        """
        recommendations = []
        seen_drugs = set()
        
        for allele in alleles:
            key = (gene, allele)
            if key in self.gene_allele_to_drugs:
                for drug_info in self.gene_allele_to_drugs[key]:
                    drug = drug_info['drug']
                    
                    # Avoid duplicate recommendations but track allele source
                    if drug in seen_drugs:
                        continue
                    seen_drugs.add(drug)
                    
                    # Categorize recommendation
                    direction = drug_info.get('direction', '').lower()
                    recommendation = {
                        'drug': drug,
                        'direction': direction,
                        'allele': allele,  # Track which allele this comes from
                        'reason': drug_info.get('pk_pd', 'Altered drug metabolism'),
                        'pmid': drug_info.get('pmid', ''),
                        'notes': drug_info.get('notes', '')
                    }
                    recommendations.append(recommendation)
        
        return recommendations
    
    def get_side_effects(self, gene: str, alleles: List[str]) -> List[Dict[str, Any]]:
        """Get potential side effects based on alleles
        
        Args:
            gene: Gene name
            alleles: List of alleles
        
        Returns:
            List of side effect warnings
        """
        side_effects = []
        seen = set()
        
        for allele in alleles:
            key = (gene, allele)
            if key in self.gene_allele_to_phenotypes:
                for info in self.gene_allele_to_phenotypes[key]:
                    se = info.get('side_effect', '')
                    if se and se not in seen:
                        seen.add(se)
                        side_effects.append({
                            'effect': se,
                            'phenotype': info.get('phenotype', ''),
                            'allele': allele,
                            'drug': info.get('drug', '')
                        })
        
        return side_effects
    
    def get_affected_drugs(self, gene: str) -> List[str]:
        """Get list of drugs affected by this gene
        
        Args:
            gene: Gene name
        
        Returns:
            List of drug names
        """
        return sorted(list(self.gene_to_drugs.get(gene, set())))
    
    def _direction_to_action(self, direction: str) -> str:
        """Convert direction string to actionable recommendation
        
        Args:
            direction: Direction of effect (e.g., 'decreased metabolism')
        
        Returns:
            Action recommendation
        """
        direction_lower = direction.lower()
        
        if any(term in direction_lower for term in ['increased risk', 'increased adverse', 'increased side', 'increased toxicity']):
            return "AVOID"
        elif any(term in direction_lower for term in ['decreased', 'reduced', 'loss', 'poor']):
            return "ADJUST DOSE"
        elif any(term in direction_lower for term in ['increased metabolism', 'increased efficacy', 'normal']):
            return "NO ADJUSTMENT"
        else:
            return "CONSULT SPECIALIST"
    
    def get_clinical_evidence(self, gene: str, alleles: List[str], diplotype: str = "") -> List[Dict[str, Any]]:
        """Get detailed clinical evidence and notes for alleles
        
        Args:
            gene: Gene name
            alleles: List of alleles
            diplotype: Optional diplotype string (e.g., '*1/*2')
        
        Returns:
            List of clinical evidence with detailed notes
        """
        evidence_list = []
        seen = set()
        
        # 1. Try diplotype lookup first (if provided)
        if diplotype:
            diplotype_variants = self._generate_diplotype_variants(diplotype)
            for dip_var in diplotype_variants:
                key = (gene, dip_var)
                self._add_evidence_from_key(key, dip_var, seen, evidence_list)
        
        # 2. Try individual alleles
        for allele in alleles:
            key = (gene, allele)
            self._add_evidence_from_key(key, allele, seen, evidence_list)
            
            # 3. Try rsID-based lookup if allele contains rsID
            if 'rs' in allele:
                # Extract rsID from patterns like "rs2231142 reference (G)" or "rs2231142"
                rsid = self._extract_rsid(allele)
                if rsid:
                    rsid_key = (gene, rsid)
                    if rsid_key in self.rsid_to_data:
                        rsid_data = self.rsid_to_data[rsid_key]
                        for info in rsid_data.get('phenotypes', []):
                            note_key = info.get('notes', '')[:100]
                            if note_key and note_key not in seen:
                                seen.add(note_key)
                                evidence_list.append({
                                    'type': 'phenotype',
                                    'allele': allele,
                                    'phenotype': info.get('phenotype', ''),
                                    'side_effect': info.get('side_effect', ''),
                                    'drug': info.get('drug', ''),
                                    'notes': info.get('notes', ''),
                                    'sentence': info.get('sentence', ''),
                                    'category': info.get('category', ''),
                                    'pmid': info.get('pmid', '')
                                })
                        
                        for info in rsid_data.get('drugs', []):
                            note_key = info.get('notes', '')[:100]
                            if note_key and note_key not in seen:
                                seen.add(note_key)
                                evidence_list.append({
                                    'type': 'drug',
                                    'allele': allele,
                                    'drug': info.get('drug', ''),
                                    'direction': info.get('direction', ''),
                                    'pk_pd': info.get('pk_pd', ''),
                                    'notes': info.get('notes', ''),
                                    'sentence': info.get('sentence', ''),
                                    'category': info.get('category', ''),
                                    'pmid': info.get('pmid', '')
                                })
                        
                        for info in rsid_data.get('functional', []):
                            note_key = info.get('notes', '')[:100]
                            if note_key and note_key not in seen:
                                seen.add(note_key)
                                evidence_list.append({
                                    'type': 'functional',
                                    'allele': allele,
                                    'drug': info.get('drug', ''),
                                    'functional_terms': info.get('functional_terms', ''),
                                    'direction': info.get('direction', ''),
                                    'notes': info.get('notes', ''),
                                    'sentence': info.get('sentence', ''),
                                    'category': info.get('category', ''),
                                    'pmid': info.get('pmid', '')
                                })
        
        return evidence_list[:15]  # Return top 15 most relevant evidence
    
    def _extract_rsid(self, allele_str: str) -> Optional[str]:
        """Extract rsID from allele string
        
        Args:
            allele_str: Allele string potentially containing rsID
        
        Returns:
            Extracted rsID or None
        """
        import re
        match = re.search(r'(rs\d+)', allele_str)
        return match.group(1) if match else None
    
    def _generate_diplotype_variants(self, diplotype: str) -> List[str]:
        """Generate common diplotype format variations
        
        Args:
            diplotype: Diplotype string
        
        Returns:
            List of diplotype format variants
        """
        if '/' not in diplotype:
            return [diplotype]
        
        alleles = [a.strip() for a in diplotype.split('/')]
        if len(alleles) != 2:
            return [diplotype]
        
        variants = set()
        for a1, a2 in [(alleles[0], alleles[1]), (alleles[1], alleles[0])]:
            # Add as-is
            variants.add(f"{a1}/{a2}")
            # Add with spaces
            variants.add(f"{a1} / {a2}")
            variants.add(f"{a1}/ {a2}")
            variants.add(f"{a1} /{a2}")
        
        return list(variants)
    
    def _add_evidence_from_key(self, key: tuple, allele: str, seen: set, evidence_list: list) -> None:
        """Helper to add evidence from a lookup key
        
        Args:
            key: Lookup key (gene, allele)
            allele: Allele name to store in evidence
            seen: Set of seen evidence (to avoid duplicates)
            evidence_list: List to append evidence to
        """
        # Get from phenotype annotations
        if key in self.gene_allele_to_phenotypes:
            for info in self.gene_allele_to_phenotypes[key]:
                note_key = info.get('notes', '')[:100]
                if note_key and note_key not in seen:
                    seen.add(note_key)
                    evidence_list.append({
                        'type': 'phenotype',
                        'allele': allele,
                        'phenotype': info.get('phenotype', ''),
                        'side_effect': info.get('side_effect', ''),
                        'drug': info.get('drug', ''),
                        'notes': info.get('notes', ''),
                        'sentence': info.get('sentence', ''),
                        'category': info.get('category', ''),
                        'pmid': info.get('pmid', '')
                    })
        
        # Get from drug annotations
        if key in self.gene_allele_to_drugs:
            for info in self.gene_allele_to_drugs[key]:
                note_key = info.get('notes', '')[:100]
                if note_key and note_key not in seen:
                    seen.add(note_key)
                    evidence_list.append({
                        'type': 'drug',
                        'allele': allele,
                        'drug': info.get('drug', ''),
                        'direction': info.get('direction', ''),
                        'pk_pd': info.get('pk_pd', ''),
                        'notes': info.get('notes', ''),
                        'sentence': info.get('sentence', ''),
                        'category': info.get('category', ''),
                        'pmid': info.get('pmid', '')
                    })
        
        # Get from functional annotations
        if key in self.gene_allele_to_functional:
            for info in self.gene_allele_to_functional[key]:
                note_key = info.get('notes', '')[:100]
                if note_key and note_key not in seen:
                    seen.add(note_key)
                    evidence_list.append({
                        'type': 'functional',
                        'allele': allele,
                        'drug': info.get('drug', ''),
                        'functional_terms': info.get('functional_terms', ''),
                        'direction': info.get('direction', ''),
                        'notes': info.get('notes', ''),
                        'sentence': info.get('sentence', ''),
                        'category': info.get('category', ''),
                        'pmid': info.get('pmid', '')
                    })
    
    def get_diplotype_drugs(self, gene: str, diplotype: str) -> Dict[str, Any]:
        """Get drug interactions specific to a diplotype
        
        Args:
            gene: Gene name
            diplotype: Diplotype string (e.g., '*1/*2' or '*2/*1')
        
        Returns:
            Dictionary with drugs and interactions for this specific diplotype, or None if not found
        """
        if not diplotype:
            return None
        
        # Try to find drugs with this exact diplotype or reverse
        alleles = diplotype.split('/')
        if len(alleles) != 2:
            return None
        
        allele1, allele2 = alleles
        
        # Check in gene_allele_to_drugs - but we need to look for exact diplotype matches
        # Since the data is indexed by individual alleles, we can only aggregate by alleles
        # For now, return None to signal that we should use the general approach
        # In a future enhancement, we could store diplotype-specific annotations
        return None
    
    def get_diplotype_evidence(self, gene: str, diplotype: str) -> List[Dict[str, Any]]:
        """Search for diplotype-specific evidence in PharmGKB annotations and allele functionality files
        
        Args:
            gene: Gene name
            diplotype: Diplotype string (e.g., '*1/*2', '*5/*6')
        
        Returns:
            List of evidence entries that specifically mention the diplotype combination
        """
        if not diplotype or '/' not in diplotype:
            return []
        
        # Normalize diplotype format
        alleles = [a.strip() for a in diplotype.split('/')]
        if len(alleles) != 2:
            return []
        
        # Generate possible diplotype representations
        diplotype_variants = set()
        for a1, a2 in [(alleles[0], alleles[1]), (alleles[1], alleles[0])]:
            diplotype_variants.add(f"{a1}/{a2}")
            a1_clean = a1.lstrip('*')
            a2_clean = a2.lstrip('*')
            diplotype_variants.add(f"{a1_clean}/{a2_clean}")
            diplotype_variants.add(f"*{a1_clean}/*{a2_clean}")
            # Also add with spaces
            diplotype_variants.add(f"{a1} /{a2}")
            diplotype_variants.add(f"{a1}/ {a2}")
            diplotype_variants.add(f"{a1} / {a2}")
        
        evidence_list = []
        seen_sentences = set()
        
        # 1. Search var_pheno_ann.tsv
        pheno_file = self.pharmgkb_dir / "var_pheno_ann.tsv"
        if pheno_file.exists():
            try:
                df = pd.read_csv(pheno_file, sep='\t', dtype=str, na_filter=False)
                
                for _, row in df.iterrows():
                    row_gene = str(row.get('Gene', '')).strip()
                    if row_gene != gene:
                        continue
                    
                    alleles_col = str(row.get('Alleles', '')).strip()
                    sentence = str(row.get('Sentence', '')).strip()
                    
                    if not alleles_col or not sentence:
                        continue
                    
                    diplotype_found = False
                    matched_diplotype = ""
                    
                    for variant in diplotype_variants:
                        if variant in alleles_col:
                            diplotype_found = True
                            matched_diplotype = variant
                            break
                    
                    if diplotype_found and sentence not in seen_sentences:
                        seen_sentences.add(sentence)
                        evidence_list.append({
                            'source': 'var_pheno_ann',
                            'diplotype': matched_diplotype,
                            'gene': row_gene,
                            'alleles': alleles_col,
                            'drug': str(row.get('Drug(s)', '')).strip(),
                            'phenotype': str(row.get('Phenotype', '')).strip(),
                            'side_effect': str(row.get('Side effect/efficacy/other', '')).strip(),
                            'direction': str(row.get('Direction of effect', '')).strip(),
                            'category': str(row.get('Phenotype Category', '')).strip(),
                            'significance': str(row.get('Significance', '')).strip(),
                            'sentence': sentence,
                            'notes': str(row.get('Notes', '')).strip(),
                            'pmid': str(row.get('PMID', '')).strip(),
                            'comparison_alleles': str(row.get('Comparison Allele(s) or Genotype(s)', '')).strip(),
                            'metabolizer_types': str(row.get('Metabolizer types', '')).strip()
                        })
            except Exception as e:
                print(f"⚠ Error searching var_pheno_ann.tsv: {e}")
        
        # 2. Search var_drug_ann.tsv
        drug_file = self.pharmgkb_dir / "var_drug_ann.tsv"
        if drug_file.exists():
            try:
                df = pd.read_csv(drug_file, sep='\t', dtype=str, na_filter=False)
                
                for _, row in df.iterrows():
                    row_gene = str(row.get('Gene', '')).strip()
                    if row_gene != gene:
                        continue
                    
                    alleles_col = str(row.get('Alleles', '')).strip()
                    sentence = str(row.get('Sentence', '')).strip()
                    
                    if not alleles_col or not sentence:
                        continue
                    
                    diplotype_found = False
                    matched_diplotype = ""
                    
                    for variant in diplotype_variants:
                        if variant in alleles_col:
                            diplotype_found = True
                            matched_diplotype = variant
                            break
                    
                    if diplotype_found and sentence not in seen_sentences:
                        seen_sentences.add(sentence)
                        evidence_list.append({
                            'source': 'var_drug_ann',
                            'diplotype': matched_diplotype,
                            'gene': row_gene,
                            'alleles': alleles_col,
                            'drug': str(row.get('Drug(s)', '')).strip(),
                            'phenotype': '',
                            'side_effect': '',
                            'direction': str(row.get('Direction of effect', '')).strip(),
                            'category': str(row.get('Phenotype Category', '')).strip(),
                            'significance': str(row.get('Significance', '')).strip(),
                            'sentence': sentence,
                            'notes': str(row.get('Notes', '')).strip(),
                            'pmid': str(row.get('PMID', '')).strip(),
                            'comparison_alleles': str(row.get('Comparison Allele(s) or Genotype(s)', '')).strip(),
                            'metabolizer_types': str(row.get('Metabolizer types', '')).strip()
                        })
            except Exception as e:
                print(f"⚠ Error searching var_drug_ann.tsv: {e}")
        
        # 3. Search allele functionality CSV file for this gene
        func_dir = self.pharmgkb_dir.parent / "source-allele-functionality"
        func_file = func_dir / f"{gene}_allele_functionality.csv"
        
        if func_file.exists():
            try:
                import csv
                with open(func_file, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        summary = str(row.get('Summary of Findings (Required)', '')).strip()
                        allele_name = str(row.get('Allele/cDNA/rsID', '')).strip()
                        
                        if not summary:
                            continue
                        
                        # Check if any diplotype variant is mentioned in the summary
                        diplotype_found = False
                        matched_diplotype = ""
                        
                        for variant in diplotype_variants:
                            if variant in summary:
                                diplotype_found = True
                                matched_diplotype = variant
                                break
                        
                        if diplotype_found:
                            # Extract just the relevant sentence mentioning the diplotype
                            sentences = summary.split('. ')
                            relevant_sentence = summary  # default to full summary
                            
                            for sent in sentences:
                                if any(v in sent for v in diplotype_variants):
                                    relevant_sentence = sent
                                    break
                            
                            if relevant_sentence not in seen_sentences:
                                seen_sentences.add(relevant_sentence)
                                
                                # Extract PMIDs from references
                                references = str(row.get('References (Required)', '')).strip()
                                pmids = references.split(',')[0].strip() if references else ''
                                
                                evidence_list.append({
                                    'source': 'allele_functionality',
                                    'diplotype': matched_diplotype,
                                    'gene': gene,
                                    'alleles': allele_name,
                                    'drug': '',
                                    'phenotype': str(row.get('Allele Clinical Functional Status (Required)', '')).strip(),
                                    'side_effect': '',
                                    'direction': '',
                                    'category': 'Functionality',
                                    'significance': str(row.get('Strength of Evidence (Required)', '')).strip(),
                                    'sentence': relevant_sentence,
                                    'notes': '',
                                    'pmid': pmids,
                                    'comparison_alleles': '',
                                    'metabolizer_types': ''
                                })
            except Exception as e:
                print(f"⚠ Error searching allele functionality file: {e}")
        
        return evidence_list
    
    def get_gene_summary(self, gene: str, diplotype: str, variant_rsids: List[str] = None) -> Dict[str, Any]:
        """Get comprehensive gene summary for report
        
        Args:
            gene: Gene name
            diplotype: Diplotype string (e.g., '*1/*2')
            variant_rsids: Optional list of rsIDs from variants
        
        Returns:
            Dictionary with all relevant clinical information
        """
        alleles = diplotype.split('/') if diplotype else []
        
        # Add variant rsIDs to allele list for lookup
        if variant_rsids:
            alleles = alleles + variant_rsids
        
        recommendations = self.get_drug_recommendations(gene, alleles)
        diplotype_evidence = self.get_diplotype_evidence(gene, diplotype)
        
        return {
            'gene': gene,
            'diplotype': diplotype,
            'diplotype_in_db': self.diplotype_exists_in_db(gene, diplotype),
            'affected_drugs': self.get_affected_drugs(gene),
            'drug_recommendations': recommendations,
            'side_effects': self.get_side_effects(gene, alleles),
            'clinical_evidence': self.get_clinical_evidence(gene, alleles, diplotype),
            'diplotype_evidence': diplotype_evidence,
            'num_recommendations': len(recommendations),
            'num_affected_drugs': len(self.get_affected_drugs(gene))
        }


class PharmGKBCache:
    """Caching layer for PharmGKB lookups to avoid repeated parsing"""
    
    def __init__(self, loader: PharmGKBDataLoader):
        """Initialize cache with loader
        
        Args:
            loader: PharmGKBDataLoader instance
        """
        self.loader = loader
        self.cache = {}
    
    def get_gene_summary(self, gene: str, diplotype: str, variant_rsids: List[str] = None) -> Dict[str, Any]:
        """Get cached gene summary"""
        rsids_key = '_'.join(sorted(variant_rsids)) if variant_rsids else ''
        key = f"{gene}_{diplotype}_{rsids_key}"
        
        if key not in self.cache:
            self.cache[key] = self.loader.get_gene_summary(gene, diplotype, variant_rsids)
        
        return self.cache[key]
    
    def clear(self) -> None:
        """Clear cache"""
        self.cache.clear()


def load_pharmgkb_data(pharmgkb_dir: str) -> PharmGKBCache:
    """Convenience function to load and cache PharmGKB data
    
    Args:
        pharmgkb_dir: Path to pharmgkb/ directory
    
    Returns:
        PharmGKBCache instance ready for queries
    """
    loader = PharmGKBDataLoader(pharmgkb_dir)
    return PharmGKBCache(loader)


# Example usage
if __name__ == "__main__":
    # Test loader
    import sys
    
    pharmgkb_dir = Path(__file__).parent.parent / "pharmgkb"
    
    if not pharmgkb_dir.exists():
        print(f"Error: PharmGKB directory not found at {pharmgkb_dir}")
        sys.exit(1)
    
    # Load data
    cache = load_pharmgkb_data(str(pharmgkb_dir))
    
    # Test query
    print("\n" + "="*70)
    print("TEST QUERY: CYP2C19 *2/*3")
    print("="*70)
    
    summary = cache.get_gene_summary("CYP2C19", "*2/*3")
    
    print(f"\nGene: {summary['gene']}")
    print(f"Diplotype: {summary['diplotype']}")
    print(f"Affected Drugs: {len(summary['affected_drugs'])} total")
    print(f"Drug Recommendations: {len(summary['drug_recommendations'])}")
    print(f"Potential Side Effects: {len(summary['side_effects'])}")
    
    if summary['drug_recommendations']:
        print("\nTop Drug Recommendations:")
        for rec in summary['drug_recommendations'][:5]:
            print(f"  • {rec['drug']}: {rec['action']} ({rec['direction']})")
    
    if summary['side_effects']:
        print("\nPotential Side Effects:")
        for se in summary['side_effects'][:5]:
            print(f"  • {se['effect']} (with {se['drug']})")
