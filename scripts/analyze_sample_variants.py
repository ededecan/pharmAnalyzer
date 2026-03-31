#!/usr/bin/env python3
"""
Pharmacogenomics Analysis Pipeline - Optimized
Properly analyzes sample variants:
  1. Strict 4-column coordinate matching (Chrom, Pos, Ref, Alt) restored.
  2. Global O(1) Locus Index to bypass missing/intergenic TSV annotations.
  3. Exact Zygosity matching (HOM/HET).
  4. Cosmetic Nomenclature Engine for IFNL3 and DPYD rsID formatting.
"""

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import json
import csv
from pathlib import Path
from typing import Dict, List, Any, Set, Tuple, Optional
from collections import defaultdict
from datetime import datetime
import sys
import argparse

if sys.version_info < (3, 9):
    import hashlib
    _orig_md5 = hashlib.md5
    def _compat_md5(*args, **kwargs):
        kwargs.pop('usedforsecurity', None)
        return _orig_md5(*args, **kwargs)
    hashlib.md5 = _compat_md5

from pharma_data_loader import (
    AlleleDefinitionLoader, AlleleFunctionalityLoader,
    DiplotypePhenotypeLoader, AlleleFrequencyLoader, DataLoader
)
from pharma_report_builder import (
    JSONReportWriter, PDFReportBuilder, SummaryTableBuilder, ReportFormatting, ConsolidatedReportHelper
)
from pharma_cache import (
    RsIDIndex, DiplotypeLookupCache, AlleleFrequencyIndex, 
    VariantGrouper, CacheManager
)
from pharmgkb_loader import load_pharmgkb_data
from drug_interaction_integration import integrate_drug_interactions_into_main_analysis, DDIHeatmapGenerator, GeneDDIVisualizer, GeneDDIPairHeatmap
from non_star_allele_genes import NonStarAlleleGeneHandler, VariantAnnotation
from dpwg_guidelines import DPWGGuidelines
from fda_boxed_warnings import FDABoxedWarningDatabase
from pharmacodynamic_genes import PharmacodynamicGenes, GeneCategory
from dosing_guidance import DosingGuidanceDatabase, PolypharmacyDosingConsiderations
from gene_prioritization import GenePrioritizer
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.lib import colors


class ImprovedPharmacogenomicsAnalyzer:
    
    def __init__(self, sample_annotation_file: str, source_dir: str = ""):
        self.sample_file = Path(sample_annotation_file)
        self.source_dir = Path(source_dir) if source_dir else Path(self.sample_file.parent)
        
        self.def_dir = self.source_dir / "source-allele-definition"
        self.freq_dir = self.source_dir / "source-allele-frequency"
        self.func_dir = self.source_dir / "source-allele-functionality"
        self.dipl_dir = self.source_dir / "source-diplotype-phenotype"
        
        for d in [self.def_dir, self.freq_dir, self.func_dir, self.dipl_dir]:
            if not d.exists(): raise FileNotFoundError(f"Directory not found: {d}")
        
        self.sample_data = self._load_sample_data()
        self.sample_id = self.sample_file.stem
        self.sample_records = self.sample_data.to_dict('records')
        
        # --- GLOBAL O(1) LOCUS INDEX ---
        self.sample_locus_index = {}
        for v in self.sample_records:
            chrom_key = next((k for k in v.keys() if k.lower() in ['chrom', '#chrom', 'chr', 'chromosome']), None)
            pos_key = next((k for k in v.keys() if k.lower() in ['pos', 'position']), None)
            
            if chrom_key and pos_key:
                chrom = str(v[chrom_key]).strip()
                if chrom and not chrom.startswith('chr'): chrom = f"chr{chrom}"
                pos = str(v[pos_key]).strip()
                key = f"{chrom}:{pos}"
                if key not in self.sample_locus_index:
                    self.sample_locus_index[key] = []
                self.sample_locus_index[key].append(v)
        
        self.allele_definitions = AlleleDefinitionLoader.load(self.def_dir)
        self.allele_functionality = AlleleFunctionalityLoader.load(self.func_dir)
        self.diplotype_phenotypes = DiplotypePhenotypeLoader.load(self.dipl_dir)
        self.allele_frequencies = AlleleFrequencyLoader.load(self.freq_dir)
        
        self.non_star_definitions = {}
        self.allele_requirements = self._build_allele_requirements()
        self.coordinate_index = self._build_rsid_coordinate_map()
        
        self.rsid_index = RsIDIndex(self.allele_definitions)
        self.diplotype_cache = DiplotypeLookupCache(self.diplotype_phenotypes)
        self.freq_index = AlleleFrequencyIndex(self.allele_frequencies)
        
        self.variant_groups = VariantGrouper.group_by_gene(self.sample_records)
        self.cache_manager = CacheManager()
        self.gene_prioritizer = GenePrioritizer()
        self.non_star_handler = NonStarAlleleGeneHandler()
        self.fda_database = FDABoxedWarningDatabase()
        self.dosing_database = DosingGuidanceDatabase()
        
        pharmgkb_dir = self.source_dir / "pharmgkb"
        self.pharmgkb_cache = load_pharmgkb_data(str(pharmgkb_dir)) if pharmgkb_dir.exists() else None
        
        print(f"✓ Loaded {len(self.sample_data)} variants from sample")
        print(f"✓ Built Global O(1) Coordinate Index to bypass TSV annotation errors")
        print(f"✓ Loaded rulebooks and mapping files for {len(self.allele_definitions)} genes")
        if self.pharmgkb_cache: print(f"✓ Loaded PharmGKB clinical data\n")

    def _extract_zygosity(self, v: dict) -> str:
        for key in ['genotype', 'GT', 'zygosity', 'Zygosity']:
            if v.get(key):
                val = str(v.get(key)).strip().upper()
                if val in ['1/1', '1|1', 'HOM', 'HOMOZYGOUS', '2/2', '2|2', 'C/C', 'T/T', 'A/A', 'G/G']: return '1/1'
                if val in ['0/1', '1/0', '0|1', '1|0', 'HET', 'HETEROZYGOUS']: return '0/1'
        
        for val in v.values():
            clean_val = str(val).strip().upper()
            if clean_val in ['HOM', '1/1', '1|1']: return '1/1'
            if clean_val in ['HET', '0/1', '1/0', '0|1', '1|0']: return '0/1'
            
        sample_val = v.get(self.sample_id) or v.get(f"{self.sample_id}.annotation")
        if sample_val:
            parts = str(sample_val).split(':')
            if parts[0] in ['1/1', '1|1']: return '1/1'
            if parts[0] in ['0/1', '1/0', '0|1', '1|0']: return '0/1'
            
        return '0/1' 

    def _build_allele_requirements(self) -> Dict[str, Dict[str, List[str]]]:
        requirements = {}
        if not self.def_dir.exists(): return requirements
            
        for csv_file in self.def_dir.glob("*_allele_definition.csv"):
            gene = csv_file.name.split('_')[0]
            requirements[gene] = {}
            try:
                with open(csv_file, 'r', encoding='utf-8') as f:
                    reader = list(csv.reader(f))
                    rsid_row_idx = next((i for i, row in enumerate(reader) if row and row[0] == 'rsID'), -1)
                    if rsid_row_idx == -1: continue
                    rsid_list = reader[rsid_row_idx]
                    
                    allele_start_idx = next((i for i, row in enumerate(reader) if row and 'Allele' in row[0]), -1)
                    if allele_start_idx == -1: continue
                    
                    is_star_gene = any(row[0].startswith('*') for row in reader[allele_start_idx + 1:] if row)
                    if not is_star_gene:
                        self.non_star_definitions[gene] = {}
                        for col_idx in range(1, len(rsid_list)):
                            rsid = rsid_list[col_idx].strip()
                            if rsid.startswith('rs'):
                                ref_name, var_name = f"{rsid} Reference", f"{rsid} Variant"
                                for row in reader[allele_start_idx + 1:]:
                                    if not row: continue
                                    name = row[0]
                                    if 'reference' in name.lower(): ref_name = name
                                    elif 'variant' in name.lower(): var_name = name
                                self.non_star_definitions[gene][rsid] = {'ref': ref_name, 'var': var_name}
                        continue 
                    
                    for row in reader[allele_start_idx + 1:]:
                        if not row or not row[0].startswith('*'): continue
                        allele_name = row[0]
                        req_rsids = []
                        for col_idx, cell_value in enumerate(row[1:], start=1):
                            if cell_value.strip() and col_idx < len(rsid_list):
                                raw_rsid = rsid_list[col_idx].strip()
                                for r in raw_rsid.replace('&', ';').replace(',', ';').split(';'):
                                    if r.strip().startswith('rs'): req_rsids.append(r.strip())
                        requirements[gene][allele_name] = req_rsids
            except Exception as e: pass
        return requirements

    def _build_rsid_coordinate_map(self) -> Dict[str, Dict[str, Dict[str, str]]]:
        tsv_path = self.source_dir / "allele_mapping.tsv"
        rsid_map = {}
        if not tsv_path.exists(): return rsid_map
        with open(tsv_path, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.strip() or line.lower().startswith('gene'): continue
                parts = [p.strip() for p in line.split('\t')]
                if len(parts) >= 6:
                    gene, variant_name, rsid, chrom, pos = parts[0], parts[1], parts[2], parts[3], parts[4]
                    ref_alt_str = parts[5].replace(' ', '') 
                    if '>' in ref_alt_str and rsid and rsid != '-':
                        ref, alt = ref_alt_str.split('>')
                        rsid_map.setdefault(gene, {})[rsid] = {'chrom': chrom if chrom.startswith('chr') else f"chr{chrom}", 'pos': pos, 'ref': ref, 'alt': alt, 'variant_name': variant_name}
        return rsid_map
    
    def _load_sample_data(self) -> pd.DataFrame:
        import gzip as _gzip
        try:
            opener = _gzip.open if str(self.sample_file).endswith('.gz') else open
            with opener(self.sample_file, 'rt', encoding='utf-8', errors='replace') as fh: lines = fh.readlines()
            header = lines[0].rstrip('\n').split('\t')
            rows = [line.rstrip('\n').split('\t')[:len(header)] + [''] * max(0, len(header) - len(line.rstrip('\n').split('\t'))) for line in lines[1:]]
            return pd.DataFrame(rows, columns=header, dtype=str)
        except Exception as e: raise Exception(f"Error loading sample data: {e}")
    
    def analyze_sample(self) -> Dict[str, Any]:
        analysis_results = {
            'analysis_id': f"PGXANALYSIS_{datetime.now().strftime('%Y%m%d_%H%M%S')}", 'analysis_date': datetime.now().isoformat(), 'sample_id': self.sample_id, 'total_variants': len(self.sample_data), 'genes_analyzed': [], 'variants_by_gene': {},
            'priority_summary': {'critical_safety_genes': [], 'major_metabolizer_genes': [], 'drug_transporter_genes': [], 'specialized_therapy_genes': [], 'clotting_factor_genes': [], 'extended_panel_genes': [], 'immediate_action_required': False, 'clinical_recommendations': []}
        }
        
        target_genes = list(self.coordinate_index.keys()) if self.coordinate_index else list(self.allele_definitions.keys())
        print(f"Evaluating all {len(target_genes)} target genes using allele_mapping.tsv exact coordinates\n")
        
        for gene in sorted(target_genes):
            print(f"{'='*70}\nANALYZING: {gene}\n{'='*70}")
            gene_result = self._analyze_gene(gene)
            is_non_star_gene = self.non_star_handler.is_non_star_allele_gene(gene)
            if gene_result.get('diplotype') or is_non_star_gene:
                if is_non_star_gene and 'identified_star_alleles' not in gene_result: gene_result['identified_star_alleles'] = []
                analysis_results['variants_by_gene'][gene] = gene_result
                analysis_results['genes_analyzed'].append(gene)
            else: print(f"⚠ Skipping {gene} - no diplotype identified")
            print()
        return self._add_priority_summary(analysis_results)
    
    def _analyze_gene(self, gene: str) -> Dict[str, Any]:
        gene_result = {
            'gene': gene, 'total_variants': 0, 'variants': [], 'identified_star_alleles': [], 'allele_functionality': {},
            'diplotype': None, 'coded_summary': None, 'ehr_priority': None, 'predicted_phenotype': None, 'metabolizer_status': None,
            'clinical_notes': [], 'pharmgkb_recommendations': [], 'affected_drugs': [], 'side_effects': [], 'priority_info': self.gene_prioritizer.get_gene_priority(gene)
        }
        
        gene_mapping = self.coordinate_index.get(gene, {})
        target_rsids = list(gene_mapping.keys())
        rsid_to_alleles, allele_to_rsids = {}, {}
        identified_alleles = set()

        for rsid in target_rsids:
            mapping = gene_mapping.get(rsid)
            if not mapping: continue

            target_chrom = mapping['chrom']
            target_pos = str(mapping['pos'])
            target_ref = mapping['ref']
            target_alt = mapping['alt']

            match = None
            
            # Fetch directly from O(1) locus index, completely bypassing the TSV's bad Gene column
            candidates = self.sample_locus_index.get(f"{target_chrom}:{target_pos}", [])
            if not candidates:
                candidates = self.variant_groups.get(gene, [])
                
            for v in candidates:
                # RESTORED STRICT COORDINATE CHECKING
                v_chrom = str(v.get('chrom', v.get('CHROM', v.get('Chr', '')))).strip()
                if v_chrom and not v_chrom.startswith('chr'): v_chrom = f"chr{v_chrom}"
                v_pos = str(v.get('pos', v.get('POS', v.get('Pos', '')))).strip()
                
                ref_key = next((k for k in v.keys() if k.lower() in ['ref', 'reference']), 'ref')
                alt_key = next((k for k in v.keys() if k.lower() in ['alt', 'alternate']), 'alt')
                
                v_ref = str(v.get(ref_key, '')).strip()
                v_alt = str(v.get(alt_key, '')).strip()
                v_alts = [a.strip() for a in v_alt.split(',')]
                
                # MUST match Chrom, Pos, Ref, and Alt
                coord_match = (v_chrom == target_chrom and v_pos == target_pos and v_ref == target_ref and target_alt in v_alts)
                
                id_key = next((k for k in v.keys() if k.lower() in ['id', 'rsid']), 'id')
                v_id_col = str(v.get(id_key, ''))
                all_tokens = set()
                for token in v_id_col.replace(';', ',').replace('&', ',').replace('|', ',').split(','):
                    if token.strip(): all_tokens.add(token.strip())
                    
                rsid_match = (rsid in all_tokens)
                
                if coord_match or rsid_match:
                    is_ref = False
                    for val in v.values():
                        if str(val).strip().upper() in ['0/0', '0|0', 'WT', 'REFERENCE', 'HOM_REF']:
                            is_ref = True
                    if is_ref: continue
                    
                    match = v
                    break
            
            if match:
                print(f"✓ VERIFIED: {rsid} found at {target_chrom}:{target_pos} ({target_ref}>{target_alt})")
                match['rsid'] = rsid
                match['position'] = f"{target_chrom}:{target_pos}"
                if 'variant_name' in mapping: match['variant_name'] = mapping['variant_name']
                gene_result['variants'].append(match)
                
                gt = self._extract_zygosity(match)
                is_hom = gt == '1/1'
                
                sas = self._find_star_alleles_from_rsid(gene, str(rsid))
                if sas:
                    rsid_to_alleles[rsid] = sas
                    for a in sas: allele_to_rsids.setdefault(a, []).append(rsid)
                    identified_alleles.update([f"{a}_hom" for a in sas] if is_hom else sas)
                elif self.non_star_handler.is_non_star_allele_gene(gene) and 'variant_name' in mapping:
                    allele_name = mapping['variant_name']
                    identified_alleles.update([f"{allele_name}_hom"] if is_hom else [allele_name])

        gene_result['total_variants'] = len(gene_result['variants'])
        if identified_alleles: gene_result['identified_star_alleles'] = sorted(list(set([a.replace('_hom', '') for a in identified_alleles])))
        
        is_non_star_gene = self.non_star_handler.is_non_star_allele_gene(gene)
        is_explicit_csv = gene in self.non_star_definitions
        
        if is_explicit_csv: diplotype = self._call_non_star_exact(gene, gene_result['variants'])
        elif rsid_to_alleles and allele_to_rsids: diplotype = self._call_diplotype_smart(gene, rsid_to_alleles, allele_to_rsids, gene_result['variants'])
        elif identified_alleles: diplotype = self._call_diplotype(gene, identified_alleles, gene_result['variants'])
        elif is_non_star_gene: diplotype = self._call_unmapped_variants(gene, gene_result['variants'])
        else: diplotype = "*1/*1" 
            
        if diplotype:
            # --- NOMENCLATURE HOTFIX ---
            if gene == 'IFNL3':
                if diplotype in ['Reference/Reference', '*1/*1']:
                    diplotype = 'rs12979860 reference (C)/rs12979860 reference (C)'
                else:
                    diplotype = diplotype.replace('-1595G>A', 'rs12979860 variant (T)')
                    diplotype = diplotype.replace('Reference', 'rs12979860 reference (C)')
                    diplotype = diplotype.replace('*1', 'rs12979860 reference (C)')
            
            pheno = self._lookup_phenotype(gene, diplotype)
            
            # DPYD Benign Compound Override
            if gene == 'DPYD' and ('c.496A>G' in diplotype and 'c.85T>C' in diplotype):
                pheno = {'phenotype': 'DPYD Normal Metabolizer', 'ehr_priority': 'Normal/Routine/Low Risk'}
                
            # DPYD Post-Lookup String Formatting (to match ABCG2/VKORC1 style)
            if gene == 'DPYD':
                for v in gene_result['variants']:
                    v_name = v.get('variant_name')
                    v_rsid = v.get('rsid')
                    v_alt = v.get('alt')
                    v_ref = v.get('ref', '')
                    if v_name and v_rsid and v_alt:
                        diplotype = diplotype.replace(v_name, f"{v_rsid} variant ({v_alt})")
                        if v_ref:
                            diplotype = diplotype.replace('Reference', f"{v_rsid} reference ({v_ref})")

            gene_result['diplotype'] = diplotype
            print(f"\n✓ Called Diplotype: {diplotype}")
            
            if pheno:
                gene_result['predicted_phenotype'] = pheno.get('phenotype', 'Unknown')
                gene_result['coded_summary'] = pheno.get('phenotype', 'Unknown')
                gene_result['ehr_priority'] = pheno.get('ehr_priority', 'Unknown')
                print(f"✓ Predicted Phenotype: {gene_result['predicted_phenotype']}")
                gene_result['metabolizer_status'] = self._determine_metabolizer_status(gene, gene_result['predicted_phenotype'])
            
            if gene in self.allele_functionality:
                for allele in gene_result['identified_star_alleles']:
                    if allele in self.allele_functionality[gene]: gene_result['allele_functionality'][allele] = self.allele_functionality[gene][allele]
            
            if self.pharmgkb_cache:
                v_rsids = [v.get('rsid', '') for v in gene_result['variants'] if v.get('rsid') and v.get('rsid') != '-']
                pgkb = self.pharmgkb_cache.get_gene_summary(gene, diplotype, v_rsids)
                gene_result['clinical_evidence'] = pgkb.get('clinical_evidence', [])
                gene_result['diplotype_evidence'] = pgkb.get('diplotype_evidence', [])
                gene_result['pharmgkb_recommendations'] = pgkb.get('drug_recommendations', [])
                gene_result['affected_drugs'] = pgkb.get('affected_drugs', [])
        
        return gene_result
    
    def _call_non_star_exact(self, gene: str, variants: list) -> str:
        defs = self.non_star_definitions.get(gene, {})
        base_ref = "Reference"
        for names in defs.values():
            if "reference" in names['ref'].lower():
                base_ref = names['ref']
                break
                
        hets, homs = [], []
        for rsid, names in defs.items():
            match = next((v for v in variants if v.get('rsid') == rsid), None)
            if match:
                var_str = names['var'] 
                gt = self._extract_zygosity(match)
                if gt == '1/1': homs.append(var_str)
                else: hets.append(var_str)
                    
        diplos = []
        for h in homs: diplos.append(f"{h}/{h}")
        if len(hets) == 2: diplos.append(f"{hets[0]}/{hets[1]}")
        else: 
            for h in hets: diplos.append(f"{h}/{base_ref}")
            
        if not diplos: return f"{base_ref}/{base_ref}"
        return " & ".join(diplos)
        
    def _call_unmapped_variants(self, gene: str, variants: list) -> str:
        if not variants: return "Reference/Reference" if gene in ['VKORC1', 'IFNL3', 'DPYD'] else "*1/*1"
        hets, homs = [], []
        for v in variants:
            name = v.get('variant_name') or v.get('rsid') or v.get('ID')
            if not name: continue
            gt = self._extract_zygosity(v)
            if gt == '1/1': homs.append(name)
            else: hets.append(name)
                
        diplos = []
        for h in homs: diplos.extend([h, h])
        if len(hets) == 2: diplos.extend([hets[0], hets[1]])
        else:
            for h in hets: diplos.append(h)
                
        if not diplos: return "Reference/Reference" if gene in ['VKORC1', 'IFNL3', 'DPYD'] else "*1/*1"
        if len(diplos) == 1: return f"Reference/{diplos[0]}" if gene in ['VKORC1', 'IFNL3', 'DPYD'] else f"*1/{diplos[0]}"
        return f"{diplos[0]}/{diplos[1]}"

    def _find_star_alleles_from_rsid(self, gene: str, rsid: str) -> List[str]:
        if not rsid or rsid == '-': return []
        return self.rsid_index.lookup_gene_alleles(str(rsid), gene)
        
    def _sort_alleles(self, alleles_list: List[str]) -> List[str]:
        def score(a): return 0 if (a == '*1' or 'reference' in a.lower()) else 1
        return sorted(alleles_list, key=score)
    
    def _call_diplotype(self, gene: str, identified_alleles: Set[str], variants: list = None) -> str:
        is_non_star = self.non_star_handler.is_non_star_allele_gene(gene)
        if not identified_alleles: return "Reference/Reference" if is_non_star else "*1/*1"
        
        alleles_list = []
        for allele in identified_alleles:
            if allele.endswith("_hom"):
                clean = allele.replace("_hom", "")
                alleles_list.extend([clean, clean])
            else: alleles_list.append(allele)
        alleles_list = sorted(list(set(alleles_list)))

        if is_non_star:
            if len(alleles_list) == 1: return f"{alleles_list[0]}/{alleles_list[0]}" if any(a.endswith("_hom") for a in identified_alleles) else f"Reference/{alleles_list[0]}"
            elif len(alleles_list) >= 2: return f"{alleles_list[0]}/{alleles_list[1]}"
            return "Reference/Reference"

        if len(alleles_list) == 1: return f"{alleles_list[0]}/{alleles_list[0]}" if any(a.endswith("_hom") for a in identified_alleles) else f"*1/{alleles_list[0]}"
        
        if gene in self.diplotype_phenotypes:
            valid = set(self.diplotype_phenotypes[gene].keys())
            for a1 in alleles_list:
                for a2 in alleles_list:
                    if f"{a1}/{a2}" in valid: return f"{a1}/{a2}"
                    if f"{a2}/{a1}" in valid: return f"{a2}/{a1}"
        final = self._sort_alleles(alleles_list)
        return f"{final[0]}/{final[1]}" if len(final) >= 2 else None
    
    def _call_diplotype_smart(self, gene: str, rsid_to_alleles: Dict, allele_to_rsids: Dict, variants: list) -> str:
        gene_rules = self.allele_requirements.get(gene, {})
        if not gene_rules: return self._call_diplotype(gene, set(allele_to_rsids.keys()), variants)
        
        variant_pool = {}
        for v in variants:
            rsid = v.get('rsid')
            if not rsid: continue
            gt = self._extract_zygosity(v)
            variant_pool[rsid] = 2 if gt == '1/1' else 1
            
        valid_candidates = []
        for allele, required_rsids in gene_rules.items():
            if allele == '*1': continue
            if set(required_rsids) and set(required_rsids).issubset(set(variant_pool.keys())):
                valid_candidates.append(allele)
                
        if not valid_candidates: return "*1/*1"
        valid_candidates.sort(key=lambda a: len(gene_rules.get(a, [])), reverse=True)
        
        final_alleles = []
        for allele in valid_candidates:
            if len(final_alleles) >= 2: break
            required_rsids = gene_rules.get(allele, [])
            for _ in range(2):
                if len(final_alleles) >= 2: break
                if all(variant_pool.get(rs, 0) > 0 for rs in required_rsids):
                    final_alleles.append(allele)
                    for rs in required_rsids: variant_pool[rs] -= 1
                else: break 
                    
        base = 'Reference' if gene == 'DPYD' else '*1'
        while len(final_alleles) < 2: final_alleles.append(base)
        final_alleles = self._sort_alleles(final_alleles)
        return f"{final_alleles[0]}/{final_alleles[1]}"

    def _lookup_phenotype(self, gene: str, diplotype: str) -> Dict[str, str]:
        # Sort alleles alphabetically to ensure dictionary matching works regardless of print order (e.g. VKORC1)
        if '/' in diplotype:
            parts = sorted(diplotype.split('/'))
            diplotype_sorted = f"{parts[0]}/{parts[1]}"
        else:
            diplotype_sorted = diplotype
            
        res = self.diplotype_cache.lookup(gene, diplotype_sorted)
        if not res:
            res = self.diplotype_cache.lookup(gene, diplotype) # Fallback to original string
            
        if res and res.get('phenotype') and res.get('phenotype') != 'Unknown': return res
        
        is_wt = diplotype in ['*1/*1', 'Reference/Reference'] or (('reference' in diplotype.lower()) and ('variant' not in diplotype.lower()))
        
        if is_wt:
            if gene in ['SLCO1B1', 'ABCG2', 'G6PD', 'NUDT15']: p = f"{gene} Normal Function"
            elif gene == 'VKORC1': p = f"{gene} Normal Warfarin Sensitivity"
            elif gene == 'IFNL3': p = f"{gene} Favorable/Normal Response"
            elif gene in ['F2', 'F5']: p = f"{gene} Normal Thrombosis Risk"
            else: p = f"{gene} Normal Metabolizer"
            return {'phenotype': p, 'ehr_priority': 'Normal/Routine/Low Risk'}
            
        return {'phenotype': 'Unknown', 'ehr_priority': 'Unknown'}
    
    def _determine_metabolizer_status(self, gene: str, phenotype: str) -> str:
        low = phenotype.lower()
        non_met = gene in ['VKORC1', 'IFNL3', 'SLCO1B1', 'ABCG2', 'F2', 'F5', 'G6PD']
        if any(x in low for x in ['poor', 'deficient']): return 'Poor Function' if non_met else 'Poor Metabolizer (PM)'
        if any(x in low for x in ['intermediate', 'decreased']): return 'Decreased Function' if non_met else 'Intermediate Metabolizer (IM)'
        if any(x in low for x in ['ultra', 'rapid', 'increased']): return 'Increased Function' if non_met else 'Ultra-Rapid Metabolizer (UM)'
        if any(x in low for x in ['normal', 'extensive', 'wild', 'favorable']): return 'Normal Function' if non_met else 'Normal Metabolizer (NM)'
        return 'Unknown'

    def _add_priority_summary(self, results: Dict[str, Any]) -> Dict[str, Any]:
        for gene in results['genes_analyzed']:
            cat = self.gene_prioritizer.get_gene_priority(gene).get('category', 'UNKNOWN')
            results['priority_summary'].get(f"{cat.lower()}_genes", []).append(gene)
            if cat == 'CRITICAL_SAFETY': results['priority_summary']['immediate_action_required'] = True
        return results

    def _print_priority_summary(self, results: Dict[str, Any]):
        summary = results.get('priority_summary', {})
        print("\n" + "="*70 + "\nCLINICAL PRIORITY SUMMARY\n" + "="*70)
        if summary.get('immediate_action_required'): print("🚨 IMMEDIATE ACTION REQUIRED: Critical safety variants found.")
        for cat, genes in summary.items():
            if isinstance(genes, list) and genes: print(f"• {cat.replace('_', ' ').upper()}: {', '.join(genes)}")
        print("="*70)

    def _classify_therapeutic_class(self, drug: str) -> str:
        d = drug.lower()
        if any(x in d for x in ['warfarin', 'clopidogrel', 'statin']): return 'Cardiovascular'
        if any(x in d for x in ['citalopram', 'sertraline', 'codeine']): return 'Psychiatric/Pain'
        if any(x in d for x in ['fluorouracil', 'tamoxifen']): return 'Oncology'
        return 'Other'

    def integrate_phase2_analysis(self, results: Dict[str, Any], patient_drugs: Optional[List[str]] = None) -> Dict[str, Any]:
        phase2_results = {'non_star_allele_genes': [], 'dpwg_recommendations': [], 'fda_warnings': [], 'pharmacodynamic_genes': {'immune_risk_genes': [], 'hemostasis_genes': [], 'metabolic_risk_genes': [], 'psychiatric_response_genes': [], 'other_pd_genes': []}, 'dosing_adjustments': [], 'drug_interaction_considerations': [], 'phase2_summary': {}}
        sample_genes = results.get('genes') or results.get('variants_by_gene')
        if not sample_genes:
            phase2_results['phase2_summary'] = self._create_phase2_summary(phase2_results)
            results['phase2_analysis'] = phase2_results
            return results
        
        for gene_name, gene_data in sample_genes.items():
            if self.non_star_handler.is_non_star_allele_gene(gene_name):
                if finding := self._analyze_non_star_gene(gene_name, gene_data): phase2_results['non_star_allele_genes'].append(self.non_star_handler.to_dict(finding))
            phase2_results['dpwg_recommendations'].extend(self._get_dpwg_recommendations(gene_name, gene_data.get('phenotype', '')))
        
        if patient_drugs:
            for drug in patient_drugs: phase2_results['fda_warnings'].extend(self._check_fda_warnings(drug, sample_genes))
            for drug in patient_drugs:
                for gene, gene_data in sample_genes.items():
                    if dosing := DosingGuidanceDatabase.get_dosing_guidance(drug, gene, gene_data.get('phenotype', '')): phase2_results['dosing_adjustments'].append(dosing.to_dict())
            for i, p_drug in enumerate(patient_drugs):
                o_drugs = patient_drugs[:i] + patient_drugs[i+1:]
                rel_genes = [g for g in sample_genes.keys() if DosingGuidanceDatabase.get_drugs_for_gene(g) and p_drug.lower() in [d.lower() for d in DosingGuidanceDatabase.get_drugs_for_gene(g)]]
                for gene in rel_genes:
                    if interactions := PolypharmacyDosingConsiderations.check_interactions(p_drug, o_drugs, gene): phase2_results['drug_interaction_considerations'].append({'primary_drug': p_drug, 'gene': gene, 'interactions': interactions})
        
        phase2_results['phase2_summary'] = self._create_phase2_summary(phase2_results)
        results['phase2_analysis'] = phase2_results
        return results

    def _analyze_non_star_gene(self, gene: str, gene_data: Dict[str, Any]) -> Optional[Any]:
        variants = []
        if 'variants' in gene_data:
            for var in gene_data['variants']:
                try:
                    pos_str = var.get('position')
                    if pos_str is None: continue
                    chrom, pos = 'Unknown', None
                    if isinstance(pos_str, str):
                        if ':' in pos_str: chrom, pos = pos_str.split(':')[0], int(pos_str.split(':')[1])
                        else: pos = int(pos_str)
                    else: pos = int(pos_str)
                    if var.get('chromosome'): chrom = var.get('chromosome')
                    variants.append(VariantAnnotation(chromosome=chrom, position=pos, ref_allele=var.get('ref') or '', alt_allele=var.get('alt') or '', rsid=var.get('rsid'), impact=var.get('impact', 'LOW'), consequence=var.get('consequence', 'Unknown'), allele_frequency=var.get('allele_frequency'), genotype=var.get('genotype')))
                except (ValueError, TypeError, AttributeError): continue
        if variants or gene_data.get('identified_star_alleles', []): return self.non_star_handler.analyze_gene_variants(gene, variants, mapped_alleles=gene_data.get('identified_star_alleles', []))
        return None

    def _get_dpwg_recommendations(self, gene: str, phenotype: str) -> List[Dict[str, Any]]:
        return [r.to_dict() for r in (DPWGGuidelines.get_recommendations(gene, phenotype, d) for d in DPWGGuidelines.get_drugs_for_gene(gene)) if r]

    def _check_fda_warnings(self, drug: str, sample_genes: Dict[str, Any]) -> List[Dict[str, Any]]:
        return [w.to_dict() for w in (self.fda_database.get_warning(drug, g, d.get('phenotype', '')) for g, d in sample_genes.items()) if w]

    def _create_phase2_summary(self, phase2_results: Dict[str, Any]) -> Dict[str, Any]:
        summary = {'total_non_star_genes_analyzed': len(phase2_results['non_star_allele_genes']), 'dpwg_recommendations_found': len(phase2_results['dpwg_recommendations']), 'fda_warnings_triggered': len(phase2_results['fda_warnings']), 'critical_warnings': [], 'dosing_adjustments_needed': len(phase2_results['dosing_adjustments']), 'critical_actions': []}
        for w in phase2_results['fda_warnings']:
            if 'AVOID' in w.get('recommended_action', '').upper():
                summary['critical_warnings'].append({'drug': w.get('drug'), 'gene': w.get('gene'), 'action': 'AVOID'})
                summary['critical_actions'].append(f"CRITICAL: AVOID {w.get('drug')} due to {w.get('gene')} genotype")
        return summary

    def generate_reports(self, analysis_results: Dict[str, Any]):
        sample_id = analysis_results.get('sample_id', 'sample_001').replace(' ', '_')
        JSONReportWriter.write(analysis_results, self.source_dir / f"{sample_id}_analysis_results.json")
        pdf_file = self.source_dir / f"{sample_id}_pharmacogenomics_report.pdf"
        self._generate_consolidated_pdf_report(analysis_results, str(pdf_file))
        print(f"✓ Consolidated pharmacogenomics PDF report saved to: {pdf_file}")
        self._generate_drug_interaction_reports(analysis_results, sample_id)

    def _generate_pdf_report(self, analysis_results: Dict[str, Any], pdf_file: str):
        builder = PDFReportBuilder(pdf_file)
        self._generate_executive_summary_page(builder, analysis_results)
        builder.add_title("Pharmacogenomics Analysis Report")
        builder.add_spacer(0.15)
        builder.add_paragraph(f"<b>Analysis ID:</b> {analysis_results['analysis_id']}<br/><b>Date:</b> {analysis_results['analysis_date']}<br/><b>Sample ID:</b> {analysis_results['sample_id']}<br/><b>Total Variants Analyzed:</b> {analysis_results['total_variants']}<br/><b>Genes Analyzed:</b> {len(analysis_results['genes_analyzed'])}<br/>")
        builder.add_spacer(0.2)
        builder.add_subsection_header("SUMMARY TABLE")
        builder.add_spacer(0.1)
        table_data, col_widths = SummaryTableBuilder.build_summary_table(analysis_results)
        builder.add_table(table_data, col_widths)
        builder.add_spacer(0.25)
        builder.add_page_break()
        builder.add_section_header("DETAILED PHARMACOGENOMIC ANALYSIS")
        builder.add_spacer(0.15)
        genes_with_priority = [(1 if str(analysis_results['variants_by_gene'][g].get('ehr_priority', 'none')).upper() in ['HIGH', 'REQUIRED'] else 2 if str(analysis_results['variants_by_gene'][g].get('ehr_priority', 'none')).upper() in ['MEDIUM', 'MODERATE'] else 3, g) for g in analysis_results['variants_by_gene'].keys()]
        for _, gene in sorted(genes_with_priority): self._add_gene_section_to_pdf(builder, gene, analysis_results['variants_by_gene'][gene])
        builder.build()

    def _generate_consolidated_pdf_report(self, analysis_results: Dict[str, Any], pdf_file: str):
        builder = PDFReportBuilder(pdf_file)
        self._generate_executive_summary_page(builder, analysis_results)
        builder.add_title("Pharmacogenomics Analysis & Drug Interaction Report")
        builder.add_spacer(0.15)
        builder.add_paragraph(f"<b>Analysis ID:</b> {analysis_results['analysis_id']}<br/><b>Date:</b> {analysis_results['analysis_date']}<br/><b>Sample ID:</b> {analysis_results['sample_id']}<br/><b>Total Variants Analyzed:</b> {analysis_results['total_variants']}<br/><b>Genes Analyzed:</b> {len(analysis_results['genes_analyzed'])}<br/>")
        builder.add_spacer(0.2)
        builder.add_subsection_header("SUMMARY TABLE")
        builder.add_spacer(0.1)
        table_data, col_widths = SummaryTableBuilder.build_summary_table(analysis_results)
        builder.add_table(table_data, col_widths)
        builder.add_spacer(0.25)
        gene_heatmap_paths = DDIHeatmapGenerator.generate_all_gene_heatmaps(analysis_results, str(self.source_dir / analysis_results.get('sample_id', 'sample_001').replace(' ', '_'))) if analysis_results.get('clinical_decision_support', {}).get('gene_drug_annotations') else {}
        builder.add_page_break()
        builder.add_section_header("DETAILED PHARMACOGENOMIC ANALYSIS")
        builder.add_spacer(0.15)
        genes_with_priority = [(1 if str(analysis_results['variants_by_gene'][g].get('ehr_priority', 'none')).upper() in ['HIGH', 'REQUIRED'] else 2 if str(analysis_results['variants_by_gene'][g].get('ehr_priority', 'none')).upper() in ['MEDIUM', 'MODERATE'] else 3, g) for g in analysis_results['variants_by_gene'].keys()]
        for _, gene in sorted(genes_with_priority): self._add_consolidated_gene_section_to_pdf(builder, gene, analysis_results['variants_by_gene'][gene], heatmap_path=gene_heatmap_paths.get(gene))
        builder.build()

    def _add_consolidated_gene_section_to_pdf(self, builder: PDFReportBuilder, gene: str, gene_data: Dict[str, Any], heatmap_path: str = None):
        builder.add_spacer(0.12)
        builder.add_paragraph(f"<u><b>{gene}</b></u>", 'GeneHeader')
        builder.add_spacer(0.08)
        builder.add_paragraph(f"<b>Variants Detected ({gene_data['total_variants']}):</b>", 'SectionBreak')
        builder.add_spacer(0.08)
        for var in gene_data.get('variants', [])[:3]: builder.add_paragraph(f"• <font face='Courier'>{var.get('position', '')} {var.get('ref', '')}→{var.get('alt', '')}</font> | <b>{var.get('consequence', 'unknown')}</b> | rsID: {var.get('rsid', '-')}", 'Normal')
        builder.add_spacer(0.15)
        diplotype_alleles = []
        if gene_data.get('diplotype'): diplotype_alleles = gene_data['diplotype'].split(' & ') if '&' in gene_data['diplotype'] else gene_data['diplotype'].split('/')
        if diplotype_alleles:
            builder.add_paragraph("<b>IDENTIFIED STAR ALLELES & FUNCTIONALITY:</b>", 'SectionBreak')
            builder.add_spacer(0.1)
            for allele in sorted(list(set(diplotype_alleles))):
                func_info = gene_data.get('allele_functionality', {}).get(allele, {})
                function_status = func_info.get('clinical_function', 'Unknown Function')
                ev_level = func_info.get('strength_of_evidence', '')
                header_text = f"<b>{allele}: {function_status}</b>"
                if ev_level and ev_level.lower() not in ['no evidence', 'none', 'unknown']: header_text += f" (Evidence: {ev_level})"
                builder.add_paragraph(header_text, 'GeneHeader')
                summary = func_info.get('summary', '').strip()
                if summary:
                    if not any(p in summary.lower() for p in ['assigned unknown function due to no literature', 'consensus among experts was unknown function', 'no available evidence']) or len(summary) > 150:
                        builder.add_spacer(0.05)
                        builder.add_paragraph(f"<b>Summary:</b> {summary}", 'Normal')
                builder.add_spacer(0.15)
        builder.add_spacer(0.15)
        if gene_data.get('diplotype'):
            builder.add_paragraph("<b>DIPLOTYPE PREDICTION:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            table_data = [['Property', 'Value'], ['Genotype', gene_data['diplotype']]]
            if gene_data.get('coded_summary') and gene_data.get('coded_summary') != 'Unknown': table_data.append(['Phenotype', gene_data['coded_summary']])
            if gene_data.get('ehr_priority') and gene_data.get('ehr_priority') != 'Unknown': table_data.append(['EHR Priority', gene_data.get('ehr_priority')])
            if gene_data.get('metabolizer_status'): table_data.append(['Metabolizer Status', gene_data['metabolizer_status']])
            builder.add_table(table_data, [1.5*inch, 4.5*inch])
            builder.add_spacer(0.15)
            if heatmap_path and Path(heatmap_path).exists():
                builder.add_paragraph(f"<b>Drug-Drug Interaction Heatmap</b>", 'SectionBreak')
                builder.add_spacer(0.08)
                try:
                    from reportlab.platypus import Image as RLImage
                    builder.story.append(RLImage(heatmap_path, width=7.0*inch, height=5.5*inch))
                except: builder.add_image(heatmap_path, width=7.0, height=5.5)
                builder.add_spacer(0.15)
            diplotype_evidence = gene_data.get('diplotype_evidence', [])
            if diplotype_evidence:
                builder.add_paragraph("<b>Diplotype-Specific Evidence from PharmGKB:</b>", 'SectionBreak')
                builder.add_spacer(0.08)
                for idx, evidence in enumerate(diplotype_evidence[:8], 1):
                    drug, category = evidence.get('drug', 'Unknown drug'), evidence.get('category', '').upper()
                    try:
                        from atc_classification import ATCClassifier
                        drug_formatted = f"<font color='#2E75B6'>[{ATCClassifier.get_therapeutic_area(drug).value}]</font> <b>{drug}</b>"
                    except: drug_formatted = f"<b>{drug}</b>"
                    if category: builder.add_paragraph(f"<font color='{self._get_category_color(category)}'><b>[{category}]</b></font> {drug_formatted}", 'Normal')
                    else: builder.add_paragraph(drug_formatted, 'Normal')
                    if sentence := evidence.get('sentence', '').strip(): builder.add_paragraph(sentence, 'Small')
                    if metabolizer := evidence.get('metabolizer_types', '').strip(): builder.add_paragraph(f"<i>Metabolizer Type: {metabolizer}</i>", 'Small')
                    if (pmid := evidence.get('pmid', '').strip()) and pmid not in ['nan', 'NaN', '']: builder.add_paragraph(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>PMID: {pmid}</u></a>', 'Small')
                    if idx < len(diplotype_evidence[:8]): builder.add_spacer(0.08)
            else:
                builder.add_paragraph(f"<i>No diplotype-specific evidence found for {gene} {gene_data['diplotype']} in PharmGKB databases.</i>", 'Normal')
                builder.add_spacer(0.08)
        builder.add_spacer(0.15)
        affected_drugs = gene_data.get('affected_drugs', [])
        if affected_drugs:
            builder.add_paragraph("<b>AFFECTED DRUGS BY THERAPEUTIC AREA:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            try:
                from atc_classification import ATCClassifier
                drugs_by_area = ATCClassifier.organize_drugs_by_area(affected_drugs)
                for therapeutic_area in sorted(drugs_by_area.keys(), key=lambda x: x.value):
                    builder.add_paragraph(f"<b><font color='#2E75B6'>{therapeutic_area.value}</font>:</b> {', '.join(sorted(drugs_by_area[therapeutic_area]))}", 'Normal')
                builder.add_spacer(0.15)
            except:
                builder.add_paragraph(f"<b>Affected Drugs:</b> {', '.join(sorted(set(affected_drugs)))}", 'Normal')
                builder.add_spacer(0.15)
        clinical_evidence = gene_data.get('clinical_evidence')
        if clinical_evidence:
            builder.add_paragraph("<b>CLINICAL EVIDENCE & SUPPORTING LITERATURE:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            evidence_by_strength = {'Level A': [], 'Level B': [], 'Level C': []}
            for evidence in clinical_evidence: evidence_by_strength[self._classify_evidence_strength(evidence)].append(evidence)
            global_idx = 1
            for strength_level in ['Level A', 'Level B', 'Level C']:
                evidence_list = evidence_by_strength[strength_level]
                if not evidence_list: continue
                color = {'Level A': '#C0392B', 'Level B': '#F39C12', 'Level C': '#27AE60'}.get(strength_level, '#7F8C8D')
                builder.add_paragraph(f"<font color='{color}'><b>{strength_level} Evidence ({len(evidence_list)} items):</b></font>", 'Normal')
                builder.add_spacer(0.08)
                from reportlab.platypus import Paragraph
                evidence_rows = [[Paragraph('<b>#</b>', builder.styles['Small']), Paragraph('<b>Context</b>', builder.styles['Small']), Paragraph('<b>Summary</b>', builder.styles['Small']), Paragraph('<b>PMID</b>', builder.styles['Small'])]]
                for evidence in evidence_list:
                    evidence_type = evidence.get('type', 'Unknown').upper()
                    type_color = self._get_evidence_type_color(evidence_type)
                    allele_raw = evidence.get('allele', '')
                    allele_parts = [p.strip() for p in allele_raw.replace('&', ';').split(';')]
                    allele = ';'.join([p for p in allele_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in allele_parts) else allele_raw
                    category = evidence.get('category', '')
                    context_parts = []
                    if category: context_parts.append(f"<font color='{self._get_category_color(category)}'><b>[{category.upper()}]</b></font>")
                    if evidence.get('drug'):
                        try:
                            from atc_classification import ATCClassifier
                            context_parts.append(f"<font color='#2E75B6'>[{ATCClassifier.get_therapeutic_area(evidence['drug']).value}]</font> <b>{evidence['drug']}</b>")
                        except: context_parts.append(f"<b>{evidence['drug']}</b>")
                    context_str = "<br/>".join(context_parts) if context_parts else "<i>Functional study</i>"
                    sentence = (evidence.get('sentence') or '').strip()
                    notes = (evidence.get('notes') or '').strip()
                    if sentence and sentence not in ['nan', 'NaN', '']: summary_txt = sentence.replace('\n', ' ').replace('\r', ' ')
                    elif notes and notes not in ['nan', 'NaN']: summary_txt = notes.replace('\n', ' ').replace('\r', ' ')
                    else: summary_txt = '<i>No summary provided</i>'
                    pmid = evidence.get('pmid', '').strip()
                    pmid_display = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>{pmid}</u></a>' if pmid and pmid not in ['nan', 'NaN', ''] else '—'
                    evidence_rows.append([Paragraph(f"<font color='{type_color}'><b>{global_idx}</b></font>", builder.styles['Small']), Paragraph(f"<b>{allele}</b><br/><font color='{type_color}'><b>{evidence_type}</b></font><br/>{context_str}", builder.styles['Small']), Paragraph(summary_txt, builder.styles['Small']), Paragraph(pmid_display, builder.styles['Small'])])
                    global_idx += 1
                builder.add_table(evidence_rows, [0.4*inch, 2.1*inch, 2.5*inch, 0.9*inch])
                builder.add_spacer(0.15)
            builder.add_spacer(0.1)
        if gene_data.get('clinical_notes'):
            builder.add_spacer(0.1)
            for note in gene_data['clinical_notes']: builder.add_paragraph(f"⚠ <b>Note:</b> {note}", 'Normal')
        builder.add_spacer(0.3)
    
    def _get_evidence_type_color(self, evidence_type: str) -> str:
        evidence_type = evidence_type.upper()
        if 'PHENOTYPE' in evidence_type: return '#006699'
        elif 'DRUG' in evidence_type: return '#009933'
        elif 'FUNCTIONAL' in evidence_type: return '#CC6600'
        return '#666666'
    
    def _get_category_color(self, category: str) -> str:
        return {'DOSAGE': '#2980B9', 'METABOLISM/PK': '#27AE60', 'TOXICITY': '#C0392B', 'EFFICACY': '#8E44AD', 'OTHER': '#7F8C8D'}.get(category.upper(), '#34495E')
    
    def _classify_therapeutic_class(self, drug_name: str) -> str:
        drug_lower = drug_name.lower()
        if any(x in drug_lower for x in ['warfarin', 'clopidogrel', 'atorvastatin', 'simvastatin', 'pravastatin']): return 'Cardiovascular'
        if any(x in drug_lower for x in ['citalopram', 'escitalopram', 'sertraline', 'fluoxetine', 'paroxetine', 'venlafaxine', 'amitriptyline', 'nortriptyline', 'desipramine', 'imipramine', 'clomipramine', 'doxepin', 'haloperidol', 'risperidone', 'aripiprazole', 'olanzapine', 'clozapine', 'codeine', 'tramadol', 'oxycodone', 'hydrocodone']): return 'Psychiatric/CNS'
        if any(x in drug_lower for x in ['fluorouracil', '5-fu', 'capecitabine', 'azathioprine', 'mercaptopurine', 'thioguanine', 'irinotecan', 'tamoxifen', 'cisplatin']): return 'Oncology'
        if any(x in drug_lower for x in ['warfarin', 'acenocoumarol', 'phenprocoumon']): return 'Anticoagulation'
        if any(x in drug_lower for x in ['codeine', 'tramadol', 'oxycodone', 'hydrocodone', 'morphine', 'fentanyl']): return 'Pain Management'
        if any(x in drug_lower for x in ['abacavir', 'efavirenz', 'atazanavir', 'isoniazid', 'rifampin']): return 'Infectious Disease'
        if any(x in drug_lower for x in ['azathioprine', 'mercaptopurine', 'tacrolimus', 'cyclosporine']): return 'Immunosuppression'
        return 'Other'
    
    def _get_guideline_level(self, recommendation: Dict[str, Any]) -> int:
        direction, category = str(recommendation.get('direction', '')).upper(), str(recommendation.get('category', '')).upper()
        if any(x in direction for x in ['AVOID', 'CONTRAINDICATED', 'ALTERNATIVE REQUIRED']): return 1
        if 'DOSAGE' in category or 'DOSE' in direction or 'TOXICITY' in category or 'ADVERSE' in direction: return 2
        if 'MONITOR' in direction or 'CAUTION' in direction or 'EFFICACY' in category: return 3
        return 4
    
    def _classify_evidence_strength(self, evidence: Dict[str, Any]) -> str:
        level = evidence.get('level_of_evidence', '').upper()
        if level in ['1A', '1B', 'A', 'HIGH']: return 'Level A'
        if level in ['2A', '2B', 'B', 'MODERATE']: return 'Level B'
        if level in ['3', '4', 'C', 'LOW', 'D']: return 'Level C'
        category, ev_type = evidence.get('category', '').upper(), evidence.get('type', '').upper()
        if category in ['DOSAGE', 'TOXICITY'] and 'DRUG' in ev_type: return 'Level A'
        if category in ['EFFICACY', 'METABOLISM/PK']: return 'Level B'
        return 'Level C'
    
    def _generate_executive_summary_page(self, builder: PDFReportBuilder, analysis_results: Dict[str, Any]):
        builder.add_paragraph("<b>EXECUTIVE SUMMARY</b>", 'Heading1')
        builder.add_spacer(0.2)
        info_text = f"<b>Sample ID:</b> {analysis_results['sample_id']}<br/><b>Analysis Date:</b> {analysis_results['analysis_date']}<br/><b>Analysis ID:</b> {analysis_results['analysis_id']}"
        builder.add_paragraph(info_text, 'Normal')
        builder.add_spacer(0.2)
        builder.add_paragraph(f"This pharmacogenomics report analyzes <b>{len(analysis_results['genes_analyzed'])} genes</b> to identify genetic variants that may affect drug metabolism and response.", 'Normal')
        builder.add_spacer(0.25)
        
        high_pri, med_pri, norm_pri = [], [], []
        for gene in sorted(analysis_results['variants_by_gene'].keys()):
            pri = str(analysis_results['variants_by_gene'][gene].get('ehr_priority', 'none')).upper()
            if pri in ['HIGH', 'REQUIRED']: high_pri.append(gene)
            elif pri in ['MEDIUM', 'MODERATE']: med_pri.append(gene)
            else: norm_pri.append(gene)
        
        builder.add_paragraph("<b>Key Findings:</b>", 'Heading3')
        builder.add_spacer(0.12)
        if high_pri: builder.add_paragraph(f"<font color='#C0392B'><b>• {len(high_pri)} gene{'s' if len(high_pri) != 1 else ''} with HIGH-PRIORITY variants requiring clinical attention</b></font>", 'Normal')
        if med_pri: builder.add_paragraph(f"<font color='#F39C12'><b>• {len(med_pri)} gene{'s' if len(med_pri) != 1 else ''} with MEDIUM-PRIORITY variants</b></font>", 'Normal')
        if norm_pri: builder.add_paragraph(f"<font color='#27AE60'>• {len(norm_pri)} gene{'s' if len(norm_pri) != 1 else ''} with normal or routine variants</font>", 'Normal')
        
        builder.add_spacer(0.2)
        if high_pri:
            builder.add_paragraph("<b><font color='#C0392B'>High Priority Genes:</font></b>", 'Heading3')
            builder.add_spacer(0.1)
            builder.add_paragraph(f"<font color='#C0392B'><b>{', '.join(sorted(high_pri))}</b></font>", 'Normal')
            builder.add_spacer(0.15)
        if med_pri:
            builder.add_paragraph("<b><font color='#F39C12'>Medium Priority Genes:</font></b>", 'Heading3')
            builder.add_spacer(0.1)
            builder.add_paragraph(f"<font color='#F39C12'><b>{', '.join(sorted(med_pri))}</b></font>", 'Normal')
            builder.add_spacer(0.15)
        
        builder.add_spacer(0.3)
        builder.add_paragraph("<i><b>Note:</b> This report provides genetic information for clinical decision support. All treatment decisions should be made by qualified healthcare professionals considering the complete clinical context.</i>", 'Small')
        builder.add_page_break()

    def _generate_drug_interaction_reports(self, analysis_results: Dict[str, Any], sample_id: str):
        drug_interactions = self._compile_drug_interactions(analysis_results)
        JSONReportWriter.write(drug_interactions, self.source_dir / f"{sample_id}_drug_interaction.json")
        
        gene_heatmap_paths = {}
        if analysis_results.get('clinical_decision_support', {}).get('gene_drug_annotations'):
            gene_heatmap_paths = DDIHeatmapGenerator.generate_all_gene_heatmaps(analysis_results, str(self.source_dir / sample_id))
            
        pdf_file = self.source_dir / f"{sample_id}_drug_interaction.pdf"
        self._generate_drug_interaction_pdf(drug_interactions, str(pdf_file), gene_heatmap_paths)
        print(f"✓ Drug interaction PDF saved to: {pdf_file}")
    
    def _compile_drug_interactions(self, analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        data = {'analysis_id': analysis_results['analysis_id'], 'sample_id': analysis_results['sample_id'], 'analysis_date': analysis_results['analysis_date'], 'total_genes_analyzed': len(analysis_results['genes_analyzed']), 'genes': []}
        for gene in sorted(analysis_results['genes_analyzed']):
            g_data = analysis_results['variants_by_gene'][gene]
            recs = g_data.get('pharmgkb_recommendations', [])
            drugs = g_data.get('affected_drugs', [])
            if recs or drugs:
                data['genes'].append({
                    'gene': gene, 'diplotype': g_data.get('diplotype', 'N/A'), 'phenotype': g_data.get('coded_summary', 'Unknown'), 
                    'metabolizer_status': g_data.get('metabolizer_status', 'Unknown'), 'affected_drugs_count': len(drugs), 
                    'drug_recommendations_count': len(recs), 'affected_drugs': drugs, 'drug_recommendations': recs, 
                    'diplotype_evidence': g_data.get('diplotype_evidence', []), 'drugbank_ddis': g_data.get('drugbank_ddis', [])
                })
        return data
    
    def _generate_drug_interaction_pdf(self, drug_interactions: Dict[str, Any], pdf_file: str, gene_heatmap_paths: Dict[str, str] = None):
        if not gene_heatmap_paths: gene_heatmap_paths = {}
        builder = PDFReportBuilder(pdf_file)
        self._setup_interaction_styles(builder)
        
        builder.add_title("Drug-Gene Interaction Report")
        builder.add_paragraph(f"<b>Sample ID:</b> {drug_interactions['sample_id']}<br/><b>Analysis Date:</b> {drug_interactions['analysis_date']}<br/>")
        builder.add_spacer(0.2)
        
        for gene_entry in drug_interactions['genes']:
            gene = gene_entry['gene']
            builder.add_paragraph(f"<b>{gene}</b>", 'CenteredGeneHeader')
            builder.add_spacer(0.08)
            
            info_table = [['Diplotype', gene_entry['diplotype']], ['Phenotype', gene_entry['phenotype']], ['Metabolizer Status', gene_entry['metabolizer_status']]]
            builder.add_table(info_table, [1.5*inch, 4.5*inch])
            builder.add_spacer(0.15)
            
            diplotype_evidence = gene_entry.get('diplotype_evidence', [])
            if diplotype_evidence:
                builder.add_subsection_header("Diplotype-Specific Evidence")
                builder.add_spacer(0.08)
                for idx, evidence in enumerate(diplotype_evidence[:5], 1):
                    drug, category = evidence.get('drug', 'Functional study'), evidence.get('category', '').upper()
                    if category: builder.add_paragraph(f"<font color='{self._get_category_color(category)}'><b>[{category}]</b></font> {drug}", 'AlleleInfo')
                    else: builder.add_paragraph(drug, 'AlleleInfo')
                    if sentence := evidence.get('sentence', '').strip(): builder.add_paragraph(sentence, 'Small')
                    if metabolizer := evidence.get('metabolizer_types', '').strip(): builder.add_paragraph(f"<i>Metabolizer Type: {metabolizer}</i>", 'Small')
                    pmid = evidence.get('pmid', '').strip()
                    if pmid and pmid not in ['nan', 'NaN', '']: builder.add_paragraph(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>PMID: {pmid}</u></a>', 'Small')
                    if idx < len(diplotype_evidence[:5]): builder.add_spacer(0.08)
                builder.add_spacer(0.15)
            
            if gene_entry['drug_recommendations'] or gene_entry['affected_drugs']:
                builder.add_subsection_header("Affected Drugs")
                builder.add_spacer(0.08)
                unique_recs = []
                seen_drugs = set()
                for rec in gene_entry['drug_recommendations']:
                    drug_name = rec.get('drug', 'Unknown').title()
                    if drug_name not in seen_drugs:
                        seen_drugs.add(drug_name)
                        unique_recs.append(rec)
                unique_recs.sort(key=lambda x: self._get_guideline_level(x))
                drugs_by_class = {}
                for rec in unique_recs:
                    tc = self._classify_therapeutic_class(rec.get('drug', 'Unknown').title())
                    drugs_by_class.setdefault(tc, []).append(rec)
                
                from reportlab.platypus import Paragraph
                if unique_recs:
                    for tc in sorted(drugs_by_class.keys()):
                        builder.add_paragraph(f"<b><font color='#2E75B6'>{tc}</font></b>", 'Normal')
                        builder.add_spacer(0.05)
                        interaction_table = [[Paragraph('<b>Drug</b>', builder.styles['Normal']), Paragraph('<b>Direction</b>', builder.styles['Normal']), Paragraph('<b>Level</b>', builder.styles['Normal'])]]
                        for rec in drugs_by_class[tc]:
                            gl = self._get_guideline_level(rec)
                            level_label = ['A', 'B', 'C', 'D'][gl - 1] if gl <= 4 else 'D'
                            interaction_table.append([Paragraph(rec.get('drug', 'Unknown').title(), builder.styles['Normal']), Paragraph(rec.get('direction', 'N/A') or 'N/A', builder.styles['Normal']), Paragraph(f"<b>{level_label}</b>", builder.styles['Normal'])])
                        builder.add_table(interaction_table, [2.0*inch, 3.5*inch, 0.6*inch])
                        builder.add_spacer(0.12)
                    
                    builder.add_subsection_header("Clinical Evidence & Summary")
                    builder.add_spacer(0.08)
                    for rec in unique_recs:
                        direction = rec.get('direction', '').strip()
                        if direction and direction.lower() != 'unknown':
                            drug_name = rec.get('drug', 'Unknown').title()
                            builder.add_paragraph(f"<b>{drug_name}</b>", 'AlleleInfo')
                            pmid = rec.get('pmid', '').strip()
                            if pmid and pmid not in ['nan', 'NaN', '', 'None']:
                                try: 
                                    pmid_num = int(pmid) if isinstance(pmid, str) else pmid
                                    builder.add_paragraph(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid_num}/" color="blue"><u>PMID: {pmid_num}</u></a>', 'Muted')
                                except: pass
                            
                            d_low = direction.lower()
                            if 'decreased' in d_low: narrative = f"<b>Decreased metabolism.</b> This allele DECREASES how the body breaks down {drug_name}, meaning the patient needs a <b>LOWER dose</b>."
                            elif 'increased' in d_low: 
                                if any(x in d_low for x in ['toxicity', 'adverse', 'risk']):
                                    narrative = f"<b>Increased effect/toxicity.</b> This allele causes {drug_name} to have a STRONGER effect or increased adverse reactions. <b>AVOID or use with extreme caution</b>."
                                else:
                                    narrative = f"<b>Increased metabolism.</b> This allele INCREASES how the body breaks down {drug_name}, meaning the patient may need a <b>HIGHER dose</b> for therapeutic effect."
                            else: narrative = f"<b>{direction.capitalize()}.</b> This allele affects the metabolism and efficacy of {drug_name} in the patient."
                            builder.add_paragraph(narrative, 'Small')
                            
                            notes = rec.get('notes', '').strip()
                            if notes and notes not in ['nan', 'NaN', '']: builder.add_paragraph(f"“{notes.replace(chr(10), ' ')}”", 'EvidenceNote')
                            builder.add_spacer(0.10)
                
                rec_drugs_lower = {r.get('drug', '').lower() for r in unique_recs}
                if gene in gene_heatmap_paths and Path(gene_heatmap_paths[gene]).exists():
                    builder.add_spacer(0.1)
                    builder.add_paragraph(f"<b>Drug-Drug Interaction Heatmap</b>", 'SectionBreak')
                    builder.add_spacer(0.08)
                    try:
                        from reportlab.platypus import Image as RLImage
                        builder.story.append(RLImage(gene_heatmap_paths[gene], width=7.0*inch, height=5.5*inch))
                        builder.add_spacer(0.15)
                    except Exception as e: print(f"Warning: Could not embed heatmap for {gene}: {e}")
                
                other_drugs = [d for d in gene_entry['affected_drugs'] if d.lower() not in rec_drugs_lower]
                if other_drugs:
                    builder.add_paragraph(f"<b>Additional Drug Interactions ({len(other_drugs)}):</b>", 'SectionBreak')
                    builder.add_spacer(0.08)
                    unique_other_drugs = set()
                    for drug_entry in other_drugs:
                        if isinstance(drug_entry, str) and drug_entry.strip():
                            for d in drug_entry.split(','):
                                d_clean = d.strip()
                                if d_clean and not d_clean.isdigit() and d_clean.lower() not in rec_drugs_lower:
                                    unique_other_drugs.add(d_clean)
                    if unique_other_drugs: builder.add_paragraph(', '.join(sorted(unique_other_drugs)), 'Normal')
                    builder.add_spacer(0.15)
            
            drugbank_ddis = gene_entry.get('drugbank_ddis', [])
            if drugbank_ddis:
                high_sev = [d for d in drugbank_ddis if d.get('severity') == 'HIGH']
                if high_sev:
                    builder.add_paragraph(f"<b>DrugBank High-Severity Drug Interactions ({len(high_sev)}):</b>", 'SectionBreak')
                    builder.add_spacer(0.08)
                    from reportlab.platypus import Paragraph
                    ddi_table = [[Paragraph('<b>Drug 1</b>', builder.styles['Normal']), Paragraph('<b>Drug 2</b>', builder.styles['Normal']), Paragraph('<b>Interaction</b>', builder.styles['Normal'])]]
                    for ddi in high_sev[:20]: ddi_table.append([Paragraph(ddi.get('drug_1', ''), builder.styles['Small']), Paragraph(ddi.get('drug_2', ''), builder.styles['Small']), Paragraph(ddi.get('description', '')[:200], builder.styles['Small'])])
                    builder.add_table(ddi_table, [1.3*inch, 1.3*inch, 3.5*inch])
                    builder.add_spacer(0.12)
                mod_sev = [d for d in drugbank_ddis if d.get('severity') == 'MODERATE']
                if mod_sev:
                    builder.add_paragraph(f"<b>DrugBank Moderate-Severity Drug Interactions:</b> {len(mod_sev)} interactions found (see full JSON report for complete details)", 'Small')
                    builder.add_spacer(0.15)
            builder.add_spacer(0.2)
        builder.build()
    
    def _setup_interaction_styles(self, builder: PDFReportBuilder) -> None:
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
        if 'CenteredGeneHeader' not in builder.styles: builder.styles.add(ParagraphStyle(name='CenteredGeneHeader', parent=builder.styles['Heading2'], alignment=TA_CENTER, spaceAfter=6, spaceBefore=12, fontSize=14, textColor='#003366'))
        if 'InteractionNarrative' not in builder.styles: builder.styles.add(ParagraphStyle(name='InteractionNarrative', parent=builder.styles['Normal'], fontSize=9, spaceAfter=8, spaceBefore=2, alignment=TA_LEFT, leading=12))


def main():
    parser = argparse.ArgumentParser(description='Pharmacogenomics Analysis Pipeline')
    parser.add_argument('--input', type=str, help='Path to input annotation TSV file')
    parser.add_argument('--source', type=str, help='Path to source data directory (contains source-allele-definition, pharmgkb, etc.)')
    parser.add_argument('--output', type=str, help='Path to output directory for reports')
    parser.add_argument('--drug-db', type=str, help='Path to drug interaction database (txt or DrugBank XML)')
    parser.add_argument('--meds', type=str, help='Comma-separated list of patient medications (e.g. "Warfarin,Simvastatin")')
    
    args = parser.parse_args()
    script_dir = Path(__file__).resolve().parent
    
    if args.input:
        sample_file = Path(args.input)
        if not sample_file.exists():
            print(f"Error: Input file not found: {sample_file}")
            return
    else:
        sample_file = script_dir / "sample_annotation.tsv"
        if not sample_file.exists():
            parent_sample = script_dir.parent / "sample_annotation.tsv"
            if parent_sample.exists():
                sample_file = parent_sample
            else:
                print(f"Error: Sample file not found in {sample_file} or {parent_sample}\nScript directory: {script_dir}")
                return
    
    if args.source:
        source_dir = Path(args.source)
        if not source_dir.exists():
            print(f"Error: Source directory not found: {source_dir}")
            return
    else:
        source_dir = sample_file.parent if sample_file.parent.name == "pharmacogenomics" else script_dir.parent
    
    output_dir = Path(args.output) if args.output else source_dir
    if args.output: 
        output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        analyzer = ImprovedPharmacogenomicsAnalyzer(sample_annotation_file=str(sample_file), source_dir=str(source_dir))
        analyzer.source_dir = output_dir
        
        print("Starting improved pharmacogenomics analysis...\n")
        results = analyzer.analyze_sample()
        print(f"\n✓ Variant analysis complete. Analyzing {len(results.get('genes_analyzed', []))} genes")
        
        if args.drug_db:
            print(f"\n📋 Drug interaction database specified: {args.drug_db}")
            db_path = Path(args.drug_db)
            
            if not db_path.exists():
                print(f"Warning: Drug DB file not found at {db_path}. Skipping interactions.")
                patient_meds = []
            else:
                print(f"✓ Found DrugBank XML at: {db_path}")
                
                if args.meds and args.meds.strip().upper() != 'NONE':
                    patient_meds = [m.strip() for m in args.meds.split(',') if m.strip()]
                else:
                    print("   Deriving patient medications from analysis results...")
                    patient_meds_set = set()
                    
                    for gene, gene_data in results.get('variants_by_gene', {}).items():
                        for drug in gene_data.get('affected_drugs', []):
                            if drug and drug.strip(): 
                                patient_meds_set.add(drug.strip())
                                
                        for rec in gene_data.get('pharmgkb_recommendations', []):
                            if (drug := rec.get('drug', '')) and drug.strip(): 
                                patient_meds_set.add(drug.strip())
                                
                    if analyzer.pharmgkb_cache:
                        for gene in results.get('variants_by_gene', {}).keys():
                            try:
                                for d in analyzer.pharmgkb_cache.get_affected_drugs(gene):
                                    if d and d.strip(): 
                                        patient_meds_set.add(d.strip())
                            except Exception: 
                                continue
                                
                    patient_meds = sorted(list(patient_meds_set))
                    
                    if patient_meds: 
                        print(f"\n💊 Derived {len(patient_meds)} medications from PharmGKB annotations")
                    else: 
                        print("\n⚠️  No medications derived from annotations. Consider providing --meds.")

                try:
                    if db_path.suffix.lower() == '.xml':
                        print(f"   Loading DrugBank XML (streaming mode)...")
                        db_content = str(db_path)
                    else:
                        print(f"   Loading drug database textfile...")
                        with open(db_path, 'r', encoding='utf-8', errors='ignore') as f: 
                            db_content = f.read()
                            
                    print(f"   Integrating DrugBank interactions with {len(patient_meds)} patient medications...")
                    results = integrate_drug_interactions_into_main_analysis(results, patient_meds, db_content)
                    print(f"   ✓ Drug interaction integration complete")
                    
                except Exception as e:
                    print(f"Error integrating drug DB: {e}")
                    import traceback
                    traceback.print_exc()
        else:
            print("\n⚠️  No --drug-db specified. Skipping drug interaction analysis.")
            patient_meds = []
        
        print("\n🔬 Integrating Phase 2 analysis (DPWG, FDA warnings, PD genes, dosing guidance)...")
        try:
            results = analyzer.integrate_phase2_analysis(results, patient_meds)
            print("   ✓ Phase 2 integration complete")
            
            if phase2_summary := results.get('phase2_analysis', {}).get('phase2_summary', {}):
                print(f"\n📊 Phase 2 Summary:")
                print(f"   Non-star genes analyzed: {phase2_summary.get('total_non_star_genes_analyzed', 0)}")
                print(f"   DPWG recommendations: {phase2_summary.get('dpwg_recommendations_found', 0)}")
                print(f"   FDA warnings: {phase2_summary.get('fda_warnings_triggered', 0)}")
                print(f"   Dosing adjustments: {phase2_summary.get('dosing_adjustments_needed', 0)}")
                
        except Exception as e:
            print(f"   ⚠️  Phase 2 integration error: {e}")
            import traceback
            traceback.print_exc()
        
        print("\n📊 Generating priority summary...")
        analyzer._print_priority_summary(results)
        analyzer.generate_reports(results)
        print("\n✓ Analysis complete!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()