#!/usr/bin/env python3
"""
Pharmacogenomics Analysis Pipeline - Optimized
Properly analyzes sample variants:
  1. Maps variants to star alleles via rsID/consequence
  2. Calls diplotypes (allele pairs)
  3. Predicts phenotypes and metabolizer status
  4. Generates comprehensive reports
"""

# Configure matplotlib backend for non-interactive plotting (must be before pyplot import)
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import json
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

# Import optimized loaders and report builders
from pharma_data_loader import (
    AlleleDefinitionLoader, AlleleFunctionalityLoader,
    DiplotypePhenotypeLoader, AlleleFrequencyLoader
)
from pharma_report_builder import (
    JSONReportWriter, PDFReportBuilder, SummaryTableBuilder, ReportFormatting, ConsolidatedReportHelper
)
from pharma_cache import (
    RsIDIndex, DiplotypeLookupCache, AlleleFrequencyIndex, 
    VariantGrouper, CacheManager
)
from pharmgkb_loader import (
    load_pharmgkb_data
)
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
    """Properly analyzes pharmacogenomic data"""
    
    def __init__(self, 
                 sample_annotation_file: str,
                 source_dir: str = ""):
        """Initialize analyzer with sample and source data"""
        self.sample_file = Path(sample_annotation_file)
        self.source_dir = Path(source_dir) if source_dir else Path(self.sample_file.parent)
        
        # Define source directories
        self.def_dir = self.source_dir / "source-allele-definition"
        self.freq_dir = self.source_dir / "source-allele-frequency"
        self.func_dir = self.source_dir / "source-allele-functionality"
        self.dipl_dir = self.source_dir / "source-diplotype-phenotype"
        
        # Validate directories exist
        for d in [self.def_dir, self.freq_dir, self.func_dir, self.dipl_dir]:
            if not d.exists():
                raise FileNotFoundError(f"Directory not found: {d}")
        
        # Load sample data
        self.sample_data = self._load_sample_data()
        
        # Extract sample ID from filename
        self.sample_id = self.sample_file.stem
        
        # Load reference databases using specialized loaders
        self.allele_definitions = AlleleDefinitionLoader.load(self.def_dir)
        self.allele_functionality = AlleleFunctionalityLoader.load(self.func_dir)
        self.diplotype_phenotypes = DiplotypePhenotypeLoader.load(self.dipl_dir)
        self.allele_frequencies = AlleleFrequencyLoader.load(self.freq_dir)
        
        # Build optimized indices for fast lookups
        self.rsid_index = RsIDIndex(self.allele_definitions)
        self.diplotype_cache = DiplotypeLookupCache(self.diplotype_phenotypes)
        self.freq_index = AlleleFrequencyIndex(self.allele_frequencies)
        self.variant_groups = VariantGrouper.group_by_gene(self.sample_data.to_dict('records'))
        
        # Initialize cache manager
        self.cache_manager = CacheManager()
        
        # Initialize gene prioritizer for clinical importance
        self.gene_prioritizer = GenePrioritizer()
        
        # Initialize non-star-allele gene handler for VKORC1, CYP4F2, IFNL3
        self.non_star_handler = NonStarAlleleGeneHandler()
        
        # Initialize Phase 2 modules
        self.fda_database = FDABoxedWarningDatabase()
        self.dosing_database = DosingGuidanceDatabase()
        
        # Load PharmGKB data for clinical recommendations
        pharmgkb_dir = self.source_dir / "pharmgkb"
        if pharmgkb_dir.exists():
            self.pharmgkb_cache = load_pharmgkb_data(str(pharmgkb_dir))
        else:
            self.pharmgkb_cache = None
        
        print(f"✓ Loaded {len(self.sample_data)} variants from sample")
        print(f"✓ Loaded allele definitions for {len(self.allele_definitions)} genes")
        print(f"✓ Loaded allele functionality for {len(self.allele_functionality)} genes")
        print(f"✓ Loaded diplotype-phenotype mapping for {len(self.diplotype_phenotypes)} genes")
        print(f"✓ Built rsID index for {len(self.rsid_index.index)} variants")
        print(f"✓ Built diplotype cache for fast lookups")
        if self.pharmgkb_cache:
            print(f"✓ Loaded PharmGKB clinical data\n")
        else:
            print(f"⚠ PharmGKB data not available\n")
    
    def _load_sample_data(self) -> pd.DataFrame:
        """Load sample annotation TSV (supports gzipped files)"""
        try:
            # Handle both regular and gzipped files
            if str(self.sample_file).endswith('.gz'):
                import gzip
                df = pd.read_csv(self.sample_file, sep='\t', dtype=str, compression='gzip')
            else:
                df = pd.read_csv(self.sample_file, sep='\t', dtype=str)
            print(f"Loaded {len(df)} variants from {self.sample_file.name}")
            return df
        except Exception as e:
            raise Exception(f"Error loading sample data: {e}")
    
    def analyze_sample(self) -> Dict[str, Any]:
        """Analyze sample and call diplotypes (optimized with pre-grouped variants)"""
        
        analysis_results = {
            'analysis_id': f"PGXANALYSIS_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            'analysis_date': datetime.now().isoformat(),
            'sample_id': self.sample_id,
            'total_variants': len(self.sample_data),
            'genes_analyzed': [],
            'variants_by_gene': {},
            'priority_summary': {
                'critical_safety_genes': [],
                'major_metabolizer_genes': [],
                'drug_transporter_genes': [],
                'specialized_therapy_genes': [],
                'clotting_factor_genes': [],
                'extended_panel_genes': [],
                'immediate_action_required': False,
                'clinical_recommendations': []
            }
        }
        
        # Use pre-grouped variants (already grouped during __init__)
        genes_in_sample = self.variant_groups
        
        # Filter genes to only those with allele definitions
        genes_with_definitions = {gene: variants for gene, variants in genes_in_sample.items() 
                                  if gene in self.allele_definitions}
        
        skipped_genes = set(genes_in_sample.keys()) - set(genes_with_definitions.keys())
        
        print(f"Found {len(genes_in_sample)} genes in sample")
        print(f"Analyzing {len(genes_with_definitions)} genes with allele definitions")
        if skipped_genes:
            print(f"Skipping {len(skipped_genes)} genes without allele definitions: {', '.join(sorted(skipped_genes))}")
        print()
        
        # Analyze each gene (pre-grouped for optimal performance)
        for gene in sorted(genes_with_definitions.keys()):
            print(f"{'='*70}")
            print(f"ANALYZING: {gene}")
            print(f"{'='*70}")
            
            variants = genes_with_definitions[gene]
            gene_result = self._analyze_gene(gene, variants)
            
            # Include genes with diplotypes OR non-star-allele genes for Phase 2 analysis
            is_non_star_gene = self.non_star_handler.is_non_star_allele_gene(gene)
            if gene_result.get('diplotype') or is_non_star_gene:
                analysis_results['variants_by_gene'][gene] = gene_result
                analysis_results['genes_analyzed'].append(gene)
                if is_non_star_gene and not gene_result.get('diplotype'):
                    print(f"ℹ️  {gene} has no star allele mapping - will be analyzed by Phase 2 modules")
            else:
                print(f"⚠ Skipping {gene} - no diplotype identified (will not be included in report)")
            print()
        
        # Generate priority summary for analyzed genes
        analysis_results = self._add_priority_summary(analysis_results)
        
        return analysis_results
    
    def _analyze_gene(self, gene: str, variants: list) -> Dict[str, Any]:
        """Analyze variants for specific gene and call diplotype"""
        
        gene_result = {
            'gene': gene,
            'total_variants': len(variants),
            'variants': [],
            'identified_star_alleles': [],
            'allele_functionality': {},
            'diplotype': None,
            'coded_summary': None,
            'ehr_priority': None,
            'predicted_phenotype': None,
            'metabolizer_status': None,
            'clinical_notes': [],
            'pharmgkb_recommendations': [],
            'affected_drugs': [],
            'side_effects': [],
            # Add priority classification
            'priority_info': self.gene_prioritizer.get_gene_priority(gene)
        }
        
        identified_alleles = set()
        
        print(f"\nVariants ({len(variants)}):")
        
        for variant in variants:
            # Extract variant info
            chrom = variant.get('chrom', '?')
            pos = variant.get('pos', '?')
            ref = variant.get('ref', '?')
            alt = variant.get('alt', '?')
            consequence = variant.get('consequence', '-')
            impact = variant.get('impact', '-')
            rsid_raw = variant.get('existing_variation', '-')
            # Filter to only show actual rsIDs (starting with 'rs'), exclude COSV, CM, etc.
            rsid = rsid_raw if rsid_raw and rsid_raw.startswith('rs') else '-'
            clinvar = variant.get('clinvar_clnsig', '-')
            hgvsp = variant.get('hgvsp', '-')
            
            variant_info = {
                'position': f"{chrom}:{pos}",
                'ref': ref,
                'alt': alt,
                'rsid': rsid,
                'consequence': consequence,
                'impact': impact,
                'hgvsp': hgvsp,
                'clinvar': clinvar
            }
            
            print(f"  • {chrom}:{pos} {ref}→{alt}")
            print(f"    rsID: {rsid} | Consequence: {consequence} ({impact})")
            
            gene_result['variants'].append(variant_info)
            
            # Try to map to star alleles using rsID
            star_alleles = self._find_star_alleles_from_rsid(gene, rsid)
            
            if star_alleles:
                print(f"    → Mapped to star alleles: {', '.join(star_alleles)}")
                identified_alleles.update(star_alleles)
            else:
                # If no direct mapping, try to infer from consequence/impact
                inferred = self._infer_star_alleles_from_consequence(gene, consequence, impact)
                if inferred:
                    print(f"    → Inferred alleles: {', '.join(inferred)}")
                    identified_alleles.update(inferred)
                else:
                    print(f"    → No star allele mapping found")
        
        # Call diplotype from identified alleles
        if identified_alleles:
            gene_result['identified_star_alleles'] = sorted(list(identified_alleles))
            
            # Get functionality information for all identified alleles
            if gene in self.allele_functionality:
                for allele in identified_alleles:
                    if allele in self.allele_functionality[gene]:
                        gene_result['allele_functionality'][allele] = self.allele_functionality[gene][allele]
            
            # Build rsid to alleles mapping and allele to rsids mapping for smart calling
            rsid_to_alleles = {}
            allele_to_rsids = {}
            
            for variant in variants:
                rsid_raw = variant.get('existing_variation', '-')
                # Filter to only actual rsIDs (starting with 'rs')
                rsid = rsid_raw if rsid_raw and rsid_raw.startswith('rs') else '-'
                if rsid and rsid != '-':
                    star_alleles = self._find_star_alleles_from_rsid(gene, rsid)
                    if star_alleles:
                        rsid_to_alleles[rsid] = star_alleles
                        for allele in star_alleles:
                            if allele not in allele_to_rsids:
                                allele_to_rsids[allele] = []
                            allele_to_rsids[allele].append(rsid)
            
            # Use smart diplotype calling if we have rsid mapping data
            if rsid_to_alleles and allele_to_rsids:
                diplotype = self._call_diplotype_smart(gene, rsid_to_alleles, allele_to_rsids, variants)
            else:
                diplotype = self._call_diplotype(gene, identified_alleles)
            
            if diplotype:
                gene_result['diplotype'] = diplotype
                print(f"\n✓ Called Diplotype: {diplotype}")
                
                # Look up phenotype and EHR priority
                phenotype_info = self._lookup_phenotype(gene, diplotype)
                if phenotype_info:
                    gene_result['predicted_phenotype'] = phenotype_info.get('phenotype', 'Unknown')
                    gene_result['coded_summary'] = phenotype_info.get('phenotype', 'Unknown')
                    gene_result['ehr_priority'] = phenotype_info.get('ehr_priority', 'Unknown')
                    print(f"✓ Predicted Phenotype: {phenotype_info.get('phenotype', 'Unknown')}")
                    print(f"✓ EHR Priority: {phenotype_info.get('ehr_priority', 'Unknown')}")
                    
                    # Determine metabolizer status
                    metabolizer = self._determine_metabolizer_status(gene, phenotype_info.get('phenotype', 'Unknown'))
                    gene_result['metabolizer_status'] = metabolizer
                    print(f"✓ Metabolizer Status: {metabolizer}")
                
                # Get PharmGKB clinical recommendations
                if self.pharmgkb_cache:
                    # Extract rsIDs from variants for additional lookup
                    variant_rsids = [v.get('rsid', '') for v in gene_result['variants'] if v.get('rsid') and v.get('rsid') != '-']
                    pharmgkb_summary = self.pharmgkb_cache.get_gene_summary(gene, diplotype, variant_rsids)
                    # Always add clinical evidence if available
                    if pharmgkb_summary.get('clinical_evidence'):
                        gene_result['clinical_evidence'] = pharmgkb_summary.get('clinical_evidence', [])
                    # Add diplotype-specific evidence
                    if pharmgkb_summary.get('diplotype_evidence'):
                        gene_result['diplotype_evidence'] = pharmgkb_summary.get('diplotype_evidence', [])
                    # Add drug recommendations and interactions
                    if pharmgkb_summary['drug_recommendations'] or pharmgkb_summary['affected_drugs']:
                        gene_result['pharmgkb_recommendations'] = pharmgkb_summary['drug_recommendations']
                        gene_result['affected_drugs'] = pharmgkb_summary['affected_drugs']
                        gene_result['side_effects'] = pharmgkb_summary['side_effects']
                        gene_result['diplotype_in_db'] = pharmgkb_summary.get('diplotype_in_db', False)
        else:
            gene_result['clinical_notes'].append("No identifiable star alleles - manual review recommended")
            print(f"\n⚠ Warning: Could not identify star alleles for {gene}")
        
        return gene_result
    
    def _find_star_alleles_from_rsid(self, gene: str, rsid: str) -> List[str]:
        """Find star alleles by rsID lookup using indexed O(1) lookup"""
        if not rsid or rsid == '-':
            return []
        
        # Use pre-built index for O(1) lookup instead of iterating
        return self.rsid_index.lookup_gene_alleles(rsid, gene)
    
    def _infer_star_alleles_from_consequence(self, gene: str, consequence: str, impact: str) -> List[str]:
        """Infer alleles from variant consequence if no rsID mapping"""
        
        if not consequence or consequence == '-':
            return []
        
        # Get available alleles for gene
        if gene not in self.allele_functionality:
            return []
        
        available_alleles = list(self.allele_functionality[gene].keys())
        
        # Map consequence to potential alleles
        if impact == 'HIGH':
            # Loss-of-function variants likely map to * alleles with loss of function
            candidates = [a for a in available_alleles if 'Loss' in self.allele_functionality[gene][a].get('function', '')]
            return candidates if candidates else ['*2']  # Default to *2 for LOF
        
        elif impact == 'MODERATE':
            # Moderate variants likely reduced function
            candidates = [a for a in available_alleles if 'Decreased' in self.allele_functionality[gene][a].get('function', '')]
            return candidates if candidates else ['*10']  # Default *10 for moderate
        
        else:
            # Low/modifier variants - assume normal function
            candidates = [a for a in available_alleles if 'Normal' in self.allele_functionality[gene][a].get('function', '')]
            return candidates if candidates else ['*1']  # Default *1 for normal
    
    def _call_diplotype(self, gene: str, identified_alleles: Set[str]) -> str:
        """Call diplotype from identified alleles with intelligent pairing"""
        
        if not identified_alleles:
            return None
        
        alleles_list = sorted(list(identified_alleles))
        
        if len(alleles_list) == 1:
            # Only one allele identified - assume homozygous
            return f"{alleles_list[0]}/{alleles_list[0]}"
        
        # Multiple alleles identified
        # Strategy: Check diplotype-phenotype mapping for valid diplotypes
        if gene in self.diplotype_phenotypes:
            valid_diplotypes = set(self.diplotype_phenotypes[gene].keys())
            
            # Try to find valid diplotypes from identified alleles
            candidates = []
            
            # Check all pairs of alleles (including homozygous)
            for a1 in alleles_list:
                for a2 in alleles_list:
                    dipl1 = f"{a1}/{a2}"
                    dipl2 = f"{a2}/{a1}"
                    
                    if dipl1 in valid_diplotypes:
                        candidates.append((dipl1, self.diplotype_phenotypes[gene][dipl1].get('activity_score', 'N/A')))
                    elif dipl2 in valid_diplotypes:
                        candidates.append((dipl2, self.diplotype_phenotypes[gene][dipl2].get('activity_score', 'N/A')))
            
            # If valid diplotypes found, return the best one (usually first two alleles or highest activity)
            if candidates:
                # Sort by activity score if available, otherwise return first
                return candidates[0][0]
        
        # Fallback: Simple pairing of first two alleles
        if len(alleles_list) >= 2:
            return f"{alleles_list[0]}/{alleles_list[1]}"
        
        return None
    
    def _call_diplotype_smart(self, gene: str, rsid_to_alleles: Dict, allele_to_rsids: Dict, variants: list) -> str:
        """Smart diplotype calling using variant patterns and allele co-occurrence"""
        
        if not rsid_to_alleles or not allele_to_rsids:
            return None
        
        # Collect all unique alleles
        all_alleles = sorted(list(allele_to_rsids.keys()))
        
        if len(all_alleles) < 1:
            return None
        
        if len(all_alleles) == 1:
            # Only one allele - homozygous
            return f"{all_alleles[0]}/{all_alleles[0]}"
        
        # IMPROVED: Use variant frequency to select most likely alleles
        # Count how many variants support each allele
        allele_support = {}
        
        for rsid, alleles in rsid_to_alleles.items():
            # All alleles mapping to this rsID get a vote
            for allele in alleles:
                if allele in all_alleles:
                    allele_support[allele] = allele_support.get(allele, 0) + 1
        
        # Sort alleles by support (frequency)
        sorted_alleles = sorted(all_alleles, 
                               key=lambda a: allele_support.get(a, 0), 
                               reverse=True)
        
        # Prefer non-*1 alleles unless *1 is significantly better supported
        non_star1 = [a for a in sorted_alleles if a != '*1']
        has_star1 = '*1' in all_alleles
        
        # Try combinations in order of preference
        candidates = []
        
        # 1. Two non-*1 alleles with highest support
        if len(non_star1) >= 2:
            a1, a2 = non_star1[0], non_star1[1]
            candidates.append((a1, a2))
        
        # 2. One non-*1 allele with *1
        if len(non_star1) >= 1 and has_star1:
            a1 = non_star1[0]
            candidates.append((a1, '*1'))
        
        # 3. One non-*1 allele homozygous
        if len(non_star1) >= 1:
            a1 = non_star1[0]
            candidates.append((a1, a1))
        
        # 4. *1 homozygous (fallback)
        if has_star1:
            candidates.append(('*1', '*1'))
        
        # 5. Top two from all alleles
        if len(sorted_alleles) >= 2:
            candidates.append((sorted_alleles[0], sorted_alleles[1]))
        
        # Try each candidate against the database
        for a1, a2 in candidates:
            for dipl in [f"{a1}/{a2}", f"{a2}/{a1}"]:
                if gene in self.diplotype_phenotypes:
                    if dipl in self.diplotype_phenotypes[gene]:
                        return dipl
                else:
                    # No database, return the first valid combination
                    return dipl
        
        # Final fallback
        if len(sorted_alleles) >= 2:
            return f"{sorted_alleles[0]}/{sorted_alleles[1]}"
        elif len(sorted_alleles) == 1:
            return f"{sorted_alleles[0]}/{sorted_alleles[0]}"
        
        return None
    
    def _lookup_phenotype(self, gene: str, diplotype: str) -> Dict[str, str]:
        """Look up phenotype information for a diplotype using cached O(1) lookup
        
        Returns:
            {
              'phenotype': Coded Diplotype/Phenotype Summary,
              'ehr_priority': EHR Priority Notation,
            }
        """
        # Use pre-built cache for O(1) lookup with automatic reversal handling
        result = self.diplotype_cache.lookup(gene, diplotype)
        
        if result:
            return result
        
        return {
            'phenotype': 'Unknown',
            'ehr_priority': 'Unknown'
        }
    
    def _determine_metabolizer_status(self, gene: str, phenotype: str) -> str:
        """Determine metabolizer status from phenotype"""
        
        phenotype_lower = phenotype.lower()
        
        # Common phenotype patterns
        if any(term in phenotype_lower for term in ['poor', 'deficient', 'loss']):
            return 'Poor Metabolizer (PM)'
        elif any(term in phenotype_lower for term in ['intermediate', 'decreased']):
            return 'Intermediate Metabolizer (IM)'
        elif any(term in phenotype_lower for term in ['ultra', 'rapid']):
            return 'Ultra-Rapid Metabolizer (UM)'
        elif any(term in phenotype_lower for term in ['normal', 'extensive', 'wild']):
            return 'Normal Metabolizer (NM)'
        else:
            return 'Unknown'
    
    def _add_priority_summary(self, analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        """Add clinical priority classification to analysis results"""
        analyzed_genes = analysis_results['genes_analyzed']
        
        # Classify genes by priority
        for gene in analyzed_genes:
            priority_info = self.gene_prioritizer.get_gene_priority(gene)
            category = priority_info.get('category', 'UNKNOWN')
            
            # Add to appropriate priority category
            if category == 'CRITICAL_SAFETY':
                analysis_results['priority_summary']['critical_safety_genes'].append(gene)
            elif category == 'MAJOR_METABOLIZERS':
                analysis_results['priority_summary']['major_metabolizer_genes'].append(gene)
            elif category == 'DRUG_TRANSPORTERS':
                analysis_results['priority_summary']['drug_transporter_genes'].append(gene)
            elif category == 'SPECIALIZED_THERAPY':
                analysis_results['priority_summary']['specialized_therapy_genes'].append(gene)
            elif category == 'CLOTTING_FACTORS':
                analysis_results['priority_summary']['clotting_factor_genes'].append(gene)
            elif category == 'EXTENDED_PANEL':
                analysis_results['priority_summary']['extended_panel_genes'].append(gene)
        
        # Generate clinical recommendations
        critical_count = len(analysis_results['priority_summary']['critical_safety_genes'])
        major_count = len(analysis_results['priority_summary']['major_metabolizer_genes'])
        
        if critical_count > 0:
            analysis_results['priority_summary']['immediate_action_required'] = True
            analysis_results['priority_summary']['clinical_recommendations'].append({
                'priority': 'URGENT',
                'action': f"Review {critical_count} critical safety gene(s): {', '.join(analysis_results['priority_summary']['critical_safety_genes'])}",
                'rationale': 'Prevent life-threatening drug toxicity'
            })
        
        if major_count > 0:
            analysis_results['priority_summary']['clinical_recommendations'].append({
                'priority': 'HIGH',
                'action': f"Evaluate {major_count} major metabolizer gene(s): {', '.join(analysis_results['priority_summary']['major_metabolizer_genes'])}",
                'rationale': 'Optimize dosing for >40% of clinical drugs'
            })
        
        # Add summary counts
        analysis_results['priority_summary']['total_genes'] = len(analyzed_genes)
        analysis_results['priority_summary']['priority_breakdown'] = {
            'CRITICAL_SAFETY': critical_count,
            'MAJOR_METABOLIZERS': major_count,
            'DRUG_TRANSPORTERS': len(analysis_results['priority_summary']['drug_transporter_genes']),
            'SPECIALIZED_THERAPY': len(analysis_results['priority_summary']['specialized_therapy_genes']),
            'CLOTTING_FACTORS': len(analysis_results['priority_summary']['clotting_factor_genes']),
            'EXTENDED_PANEL': len(analysis_results['priority_summary']['extended_panel_genes'])
        }
        
        return analysis_results
    
    def _print_priority_summary(self, analysis_results: Dict[str, Any]):
        """Print clinical priority summary to console"""
        priority_summary = analysis_results.get('priority_summary', {})
        
        print("\n" + "="*70)
        print("CLINICAL PRIORITY SUMMARY")
        print("="*70)
        
        total_genes = priority_summary.get('total_genes', 0)
        print(f"Total genes analyzed: {total_genes}")
        
        if priority_summary.get('immediate_action_required', False):
            print("\n🚨 IMMEDIATE ACTION REQUIRED:")
            critical_genes = priority_summary.get('critical_safety_genes', [])
            print(f"   Critical safety genes found: {', '.join(critical_genes)}")
            print("   → Review for life-threatening drug toxicity prevention")
        
        breakdown = priority_summary.get('priority_breakdown', {})
        if any(breakdown.values()):
            print(f"\nPRIORITY BREAKDOWN:")
            
            if breakdown.get('CRITICAL_SAFETY', 0) > 0:
                print(f"  🔴 CRITICAL SAFETY ({breakdown['CRITICAL_SAFETY']}): {', '.join(priority_summary.get('critical_safety_genes', []))}")
            
            if breakdown.get('MAJOR_METABOLIZERS', 0) > 0:
                print(f"  🟠 MAJOR METABOLIZERS ({breakdown['MAJOR_METABOLIZERS']}): {', '.join(priority_summary.get('major_metabolizer_genes', []))}")
            
            if breakdown.get('DRUG_TRANSPORTERS', 0) > 0:
                print(f"  🟡 DRUG TRANSPORTERS ({breakdown['DRUG_TRANSPORTERS']}): {', '.join(priority_summary.get('drug_transporter_genes', []))}")
            
            if breakdown.get('SPECIALIZED_THERAPY', 0) > 0:
                print(f"  🟢 SPECIALIZED THERAPY ({breakdown['SPECIALIZED_THERAPY']}): {', '.join(priority_summary.get('specialized_therapy_genes', []))}")
            
            if breakdown.get('CLOTTING_FACTORS', 0) > 0:
                print(f"  🔵 CLOTTING FACTORS ({breakdown['CLOTTING_FACTORS']}): {', '.join(priority_summary.get('clotting_factor_genes', []))}")
            
            if breakdown.get('EXTENDED_PANEL', 0) > 0:
                print(f"  ⚪ EXTENDED PANEL ({breakdown['EXTENDED_PANEL']}): {', '.join(priority_summary.get('extended_panel_genes', []))}")
        
        recommendations = priority_summary.get('clinical_recommendations', [])
        if recommendations:
            print(f"\nCLINICAL RECOMMENDATIONS:")
            for i, rec in enumerate(recommendations, 1):
                print(f"  {i}. [{rec['priority']}] {rec['action']}")
                print(f"     Rationale: {rec['rationale']}")
        
        print("="*70)
    
    def _get_evidence_type_color(self, evidence_type: str) -> str:
        """Get color for evidence type badge
        
        Args:
            evidence_type: Type of evidence (PHENOTYPE, DRUG, FUNCTIONAL)
        
        Returns:
            Hex color code
        """
        if 'PHENOTYPE' in evidence_type:
            return '#006699'  # Blue
        elif 'DRUG' in evidence_type:
            return '#009933'  # Green
        elif 'FUNCTIONAL' in evidence_type:
            return '#CC6600'  # Orange
        else:
            return '#666666'  # Gray
    
    def integrate_phase2_analysis(self, results: Dict[str, Any], patient_drugs: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Integrate Phase 2 analysis directly into results
        - Non-star-allele genes (VKORC1, CYP4F2, IFNL3)
        - DPWG guidelines
        - FDA boxed warnings
        - Pharmacodynamic genes
        - Dosing guidance
        """
        phase2_results = {
            'non_star_allele_genes': [],
            'dpwg_recommendations': [],
            'fda_warnings': [],
            'pharmacodynamic_genes': {
                'immune_risk_genes': [],
                'hemostasis_genes': [],
                'metabolic_risk_genes': [],
                'psychiatric_response_genes': [],
                'other_pd_genes': []
            },
            'dosing_adjustments': [],
            'drug_interaction_considerations': [],
            'phase2_summary': {}
        }
        
        # Extract gene analysis (support both 'genes' and 'variants_by_gene')
        sample_genes = {}
        if 'genes' in results:
            sample_genes = results['genes']
        elif 'variants_by_gene' in results:
            sample_genes = results['variants_by_gene']
        
        if not sample_genes:
            phase2_results['phase2_summary'] = self._create_phase2_summary(phase2_results)
            results['phase2_analysis'] = phase2_results
            return results
        
        # 1. Analyze non-star-allele genes
        for gene_name, gene_data in sample_genes.items():
            if self.non_star_handler.is_non_star_allele_gene(gene_name):
                finding = self._analyze_non_star_gene(gene_name, gene_data)
                if finding:
                    phase2_results['non_star_allele_genes'].append(self.non_star_handler.to_dict(finding))
        
        # 2. Get DPWG recommendations
        for gene_name, gene_data in sample_genes.items():
            phenotype = gene_data.get('phenotype', '')
            dpwg_recs = self._get_dpwg_recommendations(gene_name, phenotype)
            phase2_results['dpwg_recommendations'].extend(dpwg_recs)
        
        # 3. Check FDA warnings
        if patient_drugs:
            for drug in patient_drugs:
                warnings = self._check_fda_warnings(drug, sample_genes)
                phase2_results['fda_warnings'].extend(warnings)
        
        # 4. Pharmacodynamic gene analysis
        for gene_name, gene_data in sample_genes.items():
            info = PharmacodynamicGenes.get_gene_info(gene_name)
            if info:
                phenotype = gene_data.get('phenotype', 'Unknown')
                entry = {
                    'gene': gene_name,
                    'phenotype': phenotype,
                    'category': info.category.value,
                    'full_name': info.full_name,
                    'biological_function': info.biological_function,
                    'clinical_relevance': info.clinical_relevance,
                    'drugs_affected': info.drug_list,
                    'phenotype_implications': info.phenotype_implications.get(phenotype, 'Not specified')
                }
                
                if info.category == GeneCategory.HYPERSENSITIVITY:
                    phase2_results['pharmacodynamic_genes']['immune_risk_genes'].append(entry)
                elif info.category == GeneCategory.HEMOSTASIS:
                    phase2_results['pharmacodynamic_genes']['hemostasis_genes'].append(entry)
                elif info.category == GeneCategory.METABOLIC:
                    phase2_results['pharmacodynamic_genes']['metabolic_risk_genes'].append(entry)
                elif info.category == GeneCategory.PSYCHIATRIC:
                    phase2_results['pharmacodynamic_genes']['psychiatric_response_genes'].append(entry)
                else:
                    phase2_results['pharmacodynamic_genes']['other_pd_genes'].append(entry)
        
        # 5. Dosing adjustments
        if patient_drugs:
            for drug in patient_drugs:
                for gene, gene_data in sample_genes.items():
                    phenotype = gene_data.get('phenotype', '')
                    dosing = DosingGuidanceDatabase.get_dosing_guidance(drug, gene, phenotype)
                    if dosing:
                        phase2_results['dosing_adjustments'].append(dosing.to_dict())
        
        # 6. Drug-drug interaction considerations
        if patient_drugs:
            for i, primary_drug in enumerate(patient_drugs):
                other_drugs = patient_drugs[:i] + patient_drugs[i+1:]
                relevant_genes = [g for g in sample_genes.keys() 
                                if DosingGuidanceDatabase.get_drugs_for_gene(g) and 
                                primary_drug.lower() in [d.lower() for d in DosingGuidanceDatabase.get_drugs_for_gene(g)]]
                
                for gene in relevant_genes:
                    interactions = PolypharmacyDosingConsiderations.check_interactions(
                        primary_drug, other_drugs, gene
                    )
                    if interactions:
                        phase2_results['drug_interaction_considerations'].append({
                            'primary_drug': primary_drug,
                            'gene': gene,
                            'interactions': interactions
                        })
        
        # 7. Create Phase 2 summary
        phase2_results['phase2_summary'] = self._create_phase2_summary(phase2_results)
        
        # Add to results
        results['phase2_analysis'] = phase2_results
        return results
    
    def _analyze_non_star_gene(self, gene: str, gene_data: Dict[str, Any]) -> Optional[Any]:
        """Analyze non-star-allele gene"""
        variants = []
        if 'variants' in gene_data:
            for var in gene_data['variants']:
                try:
                    position_str = var.get('position')
                    if position_str is None:
                        continue
                    
                    # Parse position - could be "chr16:31094032" or integer
                    chromosome = 'Unknown'
                    position = None
                    
                    if isinstance(position_str, str):
                        if ':' in position_str:
                            parts = position_str.split(':')
                            chromosome = parts[0]
                            position = int(parts[1])
                        else:
                            position = int(position_str)
                    else:
                        position = int(position_str)
                    
                    if var.get('chromosome'):
                        chromosome = var.get('chromosome')
                    
                    variant = VariantAnnotation(
                        chromosome=chromosome,
                        position=position,
                        ref_allele=var.get('ref') or '',
                        alt_allele=var.get('alt') or '',
                        rsid=var.get('rsid'),
                        impact=var.get('impact', 'LOW'),
                        consequence=var.get('consequence', 'Unknown'),
                        allele_frequency=var.get('allele_frequency'),
                        genotype=var.get('genotype')
                    )
                    variants.append(variant)
                except (ValueError, TypeError, AttributeError):
                    continue
        
        if variants:
            return self.non_star_handler.analyze_gene_variants(gene, variants)
        return None
    
    def _get_dpwg_recommendations(self, gene: str, phenotype: str) -> List[Dict[str, Any]]:
        """Get DPWG recommendations for gene-phenotype"""
        recommendations = []
        drugs = DPWGGuidelines.get_drugs_for_gene(gene)
        for drug in drugs:
            dpwg_rec = DPWGGuidelines.get_recommendations(gene, phenotype, drug)
            if dpwg_rec:
                recommendations.append(dpwg_rec.to_dict())
        return recommendations
    
    def _check_fda_warnings(self, drug: str, sample_genes: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Check FDA warnings for drug given patient genotypes"""
        warnings = []
        for gene, gene_data in sample_genes.items():
            phenotype = gene_data.get('phenotype', '')
            warning = self.fda_database.get_warning(drug, gene, phenotype)
            if warning:
                warnings.append(warning.to_dict())
        return warnings
    
    def _create_phase2_summary(self, phase2_results: Dict[str, Any]) -> Dict[str, Any]:
        """Create comprehensive Phase 2 summary"""
        summary = {
            'total_non_star_genes_analyzed': len(phase2_results['non_star_allele_genes']),
            'dpwg_recommendations_found': len(phase2_results['dpwg_recommendations']),
            'fda_warnings_triggered': len(phase2_results['fda_warnings']),
            'critical_warnings': [],
            'dosing_adjustments_needed': len(phase2_results['dosing_adjustments']),
            'pharmacodynamic_genes_detected': (
                len(phase2_results['pharmacodynamic_genes'].get('immune_risk_genes', [])) +
                len(phase2_results['pharmacodynamic_genes'].get('hemostasis_genes', [])) +
                len(phase2_results['pharmacodynamic_genes'].get('metabolic_risk_genes', [])) +
                len(phase2_results['pharmacodynamic_genes'].get('psychiatric_response_genes', []))
            ),
            'critical_actions': []
        }
        
        # Identify critical warnings
        for warning in phase2_results['fda_warnings']:
            if 'AVOID' in warning.get('recommended_action', '').upper():
                summary['critical_warnings'].append({
                    'drug': warning.get('drug'),
                    'gene': warning.get('gene'),
                    'action': 'AVOID/CONTRAINDICATED'
                })
                summary['critical_actions'].append(
                    f"CRITICAL: AVOID {warning.get('drug')} due to {warning.get('gene')} genotype"
                )
        
        # Identify dosing reductions > 50%
        for dosing in phase2_results['dosing_adjustments']:
            adjustment_pct = dosing.get('adjustment_percentage')
            if adjustment_pct is not None and adjustment_pct <= -50:
                summary['critical_actions'].append(
                    f"SIGNIFICANT DOSE REDUCTION needed for {dosing.get('drug')} "
                    f"({abs(adjustment_pct)}% reduction)"
                )
        
        return summary
    
    def generate_reports(self, analysis_results: Dict[str, Any]):
        """Generate JSON and consolidated PDF reports using optimized builders"""
        
        # Extract sample ID from results
        sample_id = analysis_results.get('sample_id', 'sample_001').replace(' ', '_')
        
        # Save JSON report with sample ID prefix
        json_file = self.source_dir / f"{sample_id}_analysis_results.json"
        JSONReportWriter.write(analysis_results, json_file)
        
        # Generate consolidated PDF report that combines analysis + drug interactions
        pdf_file = self.source_dir / f"{sample_id}_pharmacogenomics_report.pdf"
        self._generate_consolidated_pdf_report(analysis_results, str(pdf_file))
        print(f"✓ Consolidated pharmacogenomics PDF report saved to: {pdf_file}")
        
        # Save legacy JSON files for backwards compatibility
        # Generate drug interaction reports JSON
        drug_interactions = self._compile_drug_interactions(analysis_results)
        ddi_json_file = self.source_dir / f"{sample_id}_drug_interaction.json"
        JSONReportWriter.write(drug_interactions, ddi_json_file)
    
    def _generate_pdf_report(self, analysis_results: Dict[str, Any], pdf_file: str):
        """Generate comprehensive PDF report using builder"""
        
        builder = PDFReportBuilder(pdf_file)
        
        # Generate executive summary page first
        self._generate_executive_summary_page(builder, analysis_results)
        
        # Title page for detailed report
        builder.add_title("Pharmacogenomics Analysis Report")
        builder.add_spacer(0.15)
        
        # Summary metadata
        summary_text = f"""
        <b>Analysis ID:</b> {analysis_results['analysis_id']}<br/>
        <b>Date:</b> {analysis_results['analysis_date']}<br/>
        <b>Sample ID:</b> {analysis_results['sample_id']}<br/>
        <b>Total Variants Analyzed:</b> {analysis_results['total_variants']}<br/>
        <b>Genes Analyzed:</b> {len(analysis_results['genes_analyzed'])}<br/>
        """
        builder.add_paragraph(summary_text)
        builder.add_spacer(0.2)
        
        # Summary table
        builder.add_subsection_header("SUMMARY TABLE")
        builder.add_spacer(0.1)
        table_data, col_widths = SummaryTableBuilder.build_summary_table(analysis_results)
        builder.add_table(table_data, col_widths)
        builder.add_spacer(0.25)
        
        # Detailed analysis
        builder.add_page_break()
        builder.add_section_header("DETAILED PHARMACOGENOMIC ANALYSIS")
        builder.add_spacer(0.15)
        
        # Sort genes by priority (high priority first)
        genes_with_priority = []
        for gene in analysis_results['variants_by_gene'].keys():
            gene_data = analysis_results['variants_by_gene'][gene]
            priority = gene_data.get('ehr_priority', 'none')
            priority_level = 3  # default
            if priority and isinstance(priority, str):
                priority_upper = priority.upper()
                if priority_upper in ['HIGH', 'REQUIRED']:
                    priority_level = 1
                elif priority_upper in ['MEDIUM', 'MODERATE']:
                    priority_level = 2
            genes_with_priority.append((priority_level, gene))
        
        # Sort by priority level, then alphabetically
        genes_with_priority.sort(key=lambda x: (x[0], x[1]))
        
        for priority_level, gene in genes_with_priority:
            gene_data = analysis_results['variants_by_gene'][gene]
            self._add_gene_section_to_pdf(builder, gene, gene_data)
        
        builder.build()
    
    def _generate_consolidated_pdf_report(self, analysis_results: Dict[str, Any], pdf_file: str):
        """Generate comprehensive consolidated PDF report with pharma analysis + drug interactions + heatmaps"""
        
        builder = PDFReportBuilder(pdf_file)
        
        # Generate executive summary page first
        self._generate_executive_summary_page(builder, analysis_results)
        
        # Title page for detailed report
        builder.add_title("Pharmacogenomics Analysis & Drug Interaction Report")
        builder.add_spacer(0.15)
        
        # Summary metadata
        summary_text = f"""
        <b>Analysis ID:</b> {analysis_results['analysis_id']}<br/>
        <b>Date:</b> {analysis_results['analysis_date']}<br/>
        <b>Sample ID:</b> {analysis_results['sample_id']}<br/>
        <b>Total Variants Analyzed:</b> {analysis_results['total_variants']}<br/>
        <b>Genes Analyzed:</b> {len(analysis_results['genes_analyzed'])}<br/>
        """
        builder.add_paragraph(summary_text)
        builder.add_spacer(0.2)
        
        # Summary table
        builder.add_subsection_header("SUMMARY TABLE")
        builder.add_spacer(0.1)
        table_data, col_widths = SummaryTableBuilder.build_summary_table(analysis_results)
        builder.add_table(table_data, col_widths)
        builder.add_spacer(0.25)
        
        # Generate gene-specific interaction heatmaps for use in comprehensive report
        gene_heatmap_paths = {}
        cds = analysis_results.get('clinical_decision_support', {})
        if cds and cds.get('gene_drug_annotations'):
            gene_heatmap_paths = DDIHeatmapGenerator.generate_all_gene_heatmaps(
                analysis_results,
                str(self.source_dir / analysis_results.get('sample_id', 'sample_001').replace(' ', '_'))
            )
        
        # Detailed analysis -  Consolidated gene sections with heatmaps and clinical evidence
        builder.add_page_break()
        builder.add_section_header("DETAILED PHARMACOGENOMIC ANALYSIS")
        builder.add_spacer(0.15)
        
        # Sort genes by priority (high priority first)
        genes_with_priority = []
        for gene in analysis_results['variants_by_gene'].keys():
            gene_data = analysis_results['variants_by_gene'][gene]
            priority = gene_data.get('ehr_priority', 'none')
            priority_level = 3  # default
            if priority and isinstance(priority, str):
                priority_upper = priority.upper()
                if priority_upper in ['HIGH', 'REQUIRED']:
                    priority_level = 1
                elif priority_upper in ['MEDIUM', 'MODERATE']:
                    priority_level = 2
            genes_with_priority.append((priority_level, gene))
        
        # Sort by priority level, then alphabetically
        genes_with_priority.sort(key=lambda x: (x[0], x[1]))
        
        # Add comprehensive gene sections with all information
        for priority_level, gene in genes_with_priority:
            gene_data = analysis_results['variants_by_gene'][gene]
            
            # Add complete gene analysis section
            self._add_consolidated_gene_section_to_pdf(
                builder, gene, gene_data, 
                heatmap_path=gene_heatmap_paths.get(gene)
            )
        
        builder.build()
    
    def _add_consolidated_gene_section_to_pdf(self, builder: PDFReportBuilder, gene: str, 
                                             gene_data: Dict[str, Any], heatmap_path: str = None):
        """Add comprehensive consolidated gene section with analysis, heatmap, and clinical evidence"""
        
        # Gene header
        builder.add_spacer(0.12)
        builder.add_paragraph(f"<u><b>{gene}</b></u>", 'GeneHeader')
        builder.add_spacer(0.08)
        
        # Variants section with monospace font and consequence type
        builder.add_paragraph(f"<b>Variants Detected ({gene_data['total_variants']}):</b>", 'SectionBreak')
        builder.add_spacer(0.08)
        
        for var in gene_data.get('variants', [])[:3]:
            var_notation = f"<font face='Courier'>{var['position']} {var['ref']}→{var['alt']}</font>"
            consequence = var.get('consequence', 'unknown')
            rsid = var.get('rsid', '-')
            
            var_text = f"• {var_notation} | <b>{consequence}</b> | rsID: {rsid}"
            builder.add_paragraph(var_text, 'Normal')
        
        builder.add_spacer(0.15)
        
        # Identified star alleles section
        diplotype_alleles = []
        if gene_data.get('diplotype'):
            diplotype_alleles = gene_data['diplotype'].split('/')
        
        if diplotype_alleles:
            builder.add_paragraph("<b>IDENTIFIED STAR ALLELES & FUNCTIONALITY:</b>", 'SectionBreak')
            builder.add_spacer(0.1)
            
            for allele in sorted(diplotype_alleles):
                func_info = gene_data.get('allele_functionality', {}).get(allele, {})
                function_status = func_info.get('clinical_function', 'Unknown Function')
                evidence_level = func_info.get('strength_of_evidence', '')
                
                if evidence_level and evidence_level.lower() not in ['no evidence', 'none', 'unknown']:
                    allele_header = f"<b>{allele}: {function_status}</b> (Evidence: {evidence_level})"
                else:
                    allele_header = f"<b>{allele}: {function_status}</b>"
                
                builder.add_paragraph(allele_header, 'GeneHeader')
                
                summary = func_info.get('summary', '').strip()
                if summary:
                    skip_phrases = [
                        'assigned unknown function due to no literature',
                        'consensus among experts was unknown function',
                        'no available evidence'
                    ]
                    summary_lower = summary.lower()
                    is_generic = any(phrase in summary_lower for phrase in skip_phrases)
                    
                    if not is_generic or len(summary) > 150:
                        builder.add_spacer(0.05)
                        builder.add_paragraph(f"<b>Summary:</b> {summary}", 'Normal')
                
                builder.add_spacer(0.15)
        
        builder.add_spacer(0.15)
        
        # Diplotype prediction section
        if gene_data.get('diplotype'):
            builder.add_paragraph("<b>DIPLOTYPE PREDICTION:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            
            table_data = [['Property', 'Value']]
            table_data.append(['Genotype', gene_data['diplotype']])
            
            if gene_data.get('coded_summary') and gene_data.get('coded_summary') != 'Unknown':
                table_data.append(['Phenotype', gene_data['coded_summary']])
            
            if gene_data.get('ehr_priority') and gene_data.get('ehr_priority') != 'Unknown':
                table_data.append(['EHR Priority', gene_data.get('ehr_priority')])
            
            if gene_data.get('metabolizer_status'):
                table_data.append(['Metabolizer Status', gene_data['metabolizer_status']])
            
            col_widths = [1.5*inch, 4.5*inch]
            builder.add_table(table_data, col_widths)
            
            builder.add_spacer(0.15)
            
            # Add heatmap for drug interactions if available
            if heatmap_path and Path(heatmap_path).exists():
                builder.add_paragraph(f"<b>Drug-Drug Interaction Heatmap</b>", 'SectionBreak')
                builder.add_spacer(0.08)
                builder.add_image(heatmap_path, width=7.0, height=5.5)
                builder.add_spacer(0.15)
            
            # Add diplotype-specific evidence from PharmGKB
            diplotype_evidence = gene_data.get('diplotype_evidence', [])
            
            if diplotype_evidence:
                builder.add_paragraph("<b>Diplotype-Specific Evidence from PharmGKB:</b>", 'SectionBreak')
                builder.add_spacer(0.08)
                
                for idx, evidence in enumerate(diplotype_evidence[:8], 1):
                    drug = evidence.get('drug', 'Unknown drug')
                    category = evidence.get('category', '').upper()
                    
                    # Format drug with therapeutic area
                    try:
                        from atc_classification import ATCClassifier
                        therapeutic_area = ATCClassifier.get_therapeutic_area(drug)
                        drug_formatted = f"<font color='#2E75B6'>[{therapeutic_area.value}]</font> <b>{drug}</b>"
                    except:
                        drug_formatted = f"<b>{drug}</b>"
                    
                    if category:
                        cat_color = self._get_category_color(category)
                        entry_header = f"<font color='{cat_color}'><b>[{category}]</b></font> {drug_formatted}"
                    else:
                        entry_header = drug_formatted
                    
                    builder.add_paragraph(entry_header, 'Normal')
                    
                    sentence = evidence.get('sentence', '').strip()
                    if sentence:
                        builder.add_paragraph(sentence, 'Small')
                    
                    metabolizer = evidence.get('metabolizer_types', '').strip()
                    if metabolizer:
                        builder.add_paragraph(f"<i>Metabolizer Type: {metabolizer}</i>", 'Small')
                    
                    pmid = evidence.get('pmid', '').strip()
                    if pmid and pmid not in ['nan', 'NaN', '']:
                        pmid_text = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>PMID: {pmid}</u></a>'
                        builder.add_paragraph(pmid_text, 'Small')
                    
                    if idx < len(diplotype_evidence[:8]):
                        builder.add_spacer(0.08)
            else:
                builder.add_paragraph(
                    f"<i>No diplotype-specific evidence found for {gene} {gene_data['diplotype']} in PharmGKB databases.</i>",
                    'Normal'
                )
                builder.add_spacer(0.08)
        
        builder.add_spacer(0.15)
        
        # Add affected drugs organized by therapeutic area
        affected_drugs = gene_data.get('affected_drugs', [])
        if affected_drugs:
            builder.add_paragraph("<b>AFFECTED DRUGS BY THERAPEUTIC AREA:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            
            try:
                from atc_classification import ATCClassifier
                drugs_by_area = ATCClassifier.organize_drugs_by_area(affected_drugs)
                
                # Sort therapeutic areas for consistent display
                for therapeutic_area in sorted(drugs_by_area.keys(), key=lambda x: x.value):
                    drugs = sorted(drugs_by_area[therapeutic_area])
                    area_text = f"<b><font color='#2E75B6'>{therapeutic_area.value}</font>:</b> " + ", ".join(drugs)
                    builder.add_paragraph(area_text, 'Normal')
                
                builder.add_spacer(0.15)
            except Exception as e:
                # Fallback if classification fails
                drugs_text = ", ".join(sorted(set(affected_drugs)))
                builder.add_paragraph(f"<b>Affected Drugs:</b> {drugs_text}", 'Normal')
                builder.add_spacer(0.15)
        
        # Add comprehensive clinical evidence after heatmap/diplotype evidence
        if gene_data.get('clinical_evidence'):
            builder.add_paragraph("<b>CLINICAL EVIDENCE & SUPPORTING LITERATURE:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            
            # Group evidence by strength
            evidence_by_strength = {'Level A': [], 'Level B': [], 'Level C': []}
            for evidence in gene_data['clinical_evidence']:
                strength = self._classify_evidence_strength(evidence)
                evidence_by_strength[strength].append(evidence)
            
            # Display evidence grouped by strength
            global_idx = 1
            for strength_level in ['Level A', 'Level B', 'Level C']:
                evidence_list = evidence_by_strength[strength_level]
                if not evidence_list:
                    continue
                
                # Strength level header with color coding
                strength_colors = {
                    'Level A': '#C0392B',
                    'Level B': '#F39C12', 
                    'Level C': '#27AE60'
                }
                color = strength_colors.get(strength_level, '#7F8C8D')
                builder.add_paragraph(f"<font color='{color}'><b>{strength_level} Evidence ({len(evidence_list)} items):</b></font>", 'Normal')
                builder.add_spacer(0.08)
                
                # Build evidence table for this strength level
                from reportlab.platypus import Paragraph
                evidence_rows = [
                    [
                        Paragraph('<b>#</b>', builder.styles['Small']),
                        Paragraph('<b>Context</b>', builder.styles['Small']),
                        Paragraph('<b>Summary</b>', builder.styles['Small']),
                        Paragraph('<b>PMID</b>', builder.styles['Small'])
                    ]
                ]

                for evidence in evidence_list:
                    evidence_type = evidence.get('type', 'Unknown').upper()
                    type_color = self._get_evidence_type_color(evidence_type)
                    allele_raw = evidence.get('allele', '')
                    allele_parts = [p.strip() for p in allele_raw.replace('&', ';').split(';')]
                    allele = ';'.join([p for p in allele_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in allele_parts) else allele_raw
                    category = evidence.get('category', '')

                    # Build context with category badge and drug (with therapeutic area)
                    context_parts = []
                    
                    if category:
                        cat_color = self._get_category_color(category)
                        context_parts.append(f"<font color='{cat_color}'><b>[{category.upper()}]</b></font>")
                    
                    if evidence.get('drug'):
                        # Add therapeutic area classification
                        try:
                            from atc_classification import ATCClassifier
                            therapeutic_area = ATCClassifier.get_therapeutic_area(evidence['drug'])
                            context_parts.append(f"<font color='#2E75B6'>[{therapeutic_area.value}]</font> <b>{evidence['drug']}</b>")
                        except:
                            context_parts.append(f"<b>{evidence['drug']}</b>")
                    
                    if context_parts:
                        context_str = "<br/>".join(context_parts)
                    else:
                        context_str = "<i>Functional study</i>"

                    # Use Sentence field if available
                    sentence = (evidence.get('sentence') or '').strip()
                    notes = (evidence.get('notes') or '').strip()
                    
                    if sentence and sentence not in ['nan', 'NaN', '']:
                        summary_txt = sentence.replace('\n', ' ').replace('\r', ' ')
                    elif notes and notes not in ['nan', 'NaN']:
                        summary_txt = notes.replace('\n', ' ').replace('\r', ' ')
                    else:
                        summary_txt = '<i>No summary provided</i>'

                    # Make PMID clickable
                    pmid = evidence.get('pmid', '').strip()
                    if pmid and pmid not in ['nan', 'NaN', '']:
                        pmid_display = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>{pmid}</u></a>'
                    else:
                        pmid_display = '—'

                    evidence_rows.append([
                        Paragraph(f"<font color='{type_color}'><b>{global_idx}</b></font>", builder.styles['Small']),
                        Paragraph(f"<b>{allele}</b><br/><font color='{type_color}'><b>{evidence_type}</b></font><br/>{context_str}", builder.styles['Small']),
                        Paragraph(summary_txt, builder.styles['Small']),
                        Paragraph(pmid_display, builder.styles['Small'])
                    ])
                    global_idx += 1

                # Add evidence table for this strength level
                col_widths = [0.4*inch, 2.1*inch, 2.5*inch, 0.9*inch]
                builder.add_table(evidence_rows, col_widths)
                builder.add_spacer(0.15)

            builder.add_spacer(0.1)
        
        # Clinical notes
        if gene_data.get('clinical_notes'):
            builder.add_spacer(0.1)
            for note in gene_data['clinical_notes']:
                builder.add_paragraph(f"⚠ <b>Note:</b> {note}", 'Normal')
        
        builder.add_spacer(0.3)
    
    def _generate_drug_interaction_reports(self, analysis_results: Dict[str, Any], sample_id: str):
        """Generate comprehensive drug interaction reports (PDF and JSON)"""
        
        # Compile drug interaction data
        drug_interactions = self._compile_drug_interactions(analysis_results)
        
        # Generate JSON report
        json_file = self.source_dir / f"{sample_id}_drug_interaction.json"
        JSONReportWriter.write(drug_interactions, json_file)
        
        # Generate gene-specific interaction heatmaps
        gene_heatmap_paths = {}
        cds = analysis_results.get('clinical_decision_support', {})
        if cds and cds.get('gene_drug_annotations'):
            # Generate drug pair interaction heatmaps (primary visualization)
            gene_heatmap_paths = DDIHeatmapGenerator.generate_all_gene_heatmaps(
                analysis_results,
                str(self.source_dir / sample_id)
            )
        
        # Generate PDF report (use heatmaps as primary visualization)
        pdf_file = self.source_dir / f"{sample_id}_drug_interaction.pdf"
        self._generate_drug_interaction_pdf(drug_interactions, str(pdf_file), gene_heatmap_paths)
        print(f"✓ Drug interaction PDF saved to: {pdf_file}")
    
    def _compile_drug_interactions(self, analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        """Compile all drug interactions across all genes"""
        
        drug_interactions_data = {
            'analysis_id': analysis_results['analysis_id'],
            'sample_id': analysis_results['sample_id'],
            'analysis_date': analysis_results['analysis_date'],
            'total_genes_analyzed': len(analysis_results['genes_analyzed']),
            'genes': []
        }
        
        # Compile drug interactions for each gene
        for gene in sorted(analysis_results['genes_analyzed']):
            gene_data = analysis_results['variants_by_gene'][gene]
            
            # Get drug recommendations and affected drugs
            drug_recommendations = gene_data.get('pharmgkb_recommendations', [])
            affected_drugs = gene_data.get('affected_drugs', [])
            
            # Only include gene entry if it has drug interaction data
            if drug_recommendations or affected_drugs:
                gene_entry = {
                    'gene': gene,
                    'diplotype': gene_data.get('diplotype', 'N/A'),
                    'phenotype': gene_data.get('coded_summary', 'Unknown'),
                    'metabolizer_status': gene_data.get('metabolizer_status', 'Unknown'),
                    'affected_drugs_count': len(affected_drugs),
                    'drug_recommendations_count': len(drug_recommendations),
                    'affected_drugs': affected_drugs,
                    'drug_recommendations': drug_recommendations,
                    'diplotype_evidence': gene_data.get('diplotype_evidence', []),
                    'drugbank_ddis': gene_data.get('drugbank_ddis', [])
                }
                
                drug_interactions_data['genes'].append(gene_entry)
        
        return drug_interactions_data
    
    def _generate_drug_interaction_pdf(self, drug_interactions: Dict[str, Any], pdf_file: str, 
                                      gene_heatmap_paths: Dict[str, str] = None):
        """Generate comprehensive drug interaction PDF with drug pair heatmaps"""
        
        if gene_heatmap_paths is None:
            gene_heatmap_paths = {}
        
        builder = PDFReportBuilder(pdf_file)
        
        # Setup custom interaction styles
        self._setup_interaction_styles(builder)
        
        # Title
        builder.add_title("Drug-Gene Interaction Report")
        
        # Summary metadata
        summary_text = f"""
        <b>Sample ID:</b> {drug_interactions['sample_id']}<br/>
        <b>Analysis Date:</b> {drug_interactions['analysis_date']}<br/>
        """
        builder.add_paragraph(summary_text)
        builder.add_spacer(0.2)
        
        # Drug interaction section for each gene
        for gene_entry in drug_interactions['genes']:
            gene = gene_entry['gene']
            diplotype = gene_entry['diplotype']
            
            # Gene name as main header (centered)
            builder.add_paragraph(f"<b>{gene}</b>", 'CenteredGeneHeader')
            builder.add_spacer(0.08)
            
            # Gene info table
            gene_info_data = [
                ['Diplotype', gene_entry['diplotype']],
                ['Phenotype', gene_entry['phenotype']],
                ['Metabolizer Status', gene_entry['metabolizer_status']]
            ]
            builder.add_table(gene_info_data, [1.5*inch, 4.5*inch])
            builder.add_spacer(0.15)
            
            # Add diplotype-specific evidence if available
            diplotype_evidence = gene_entry.get('diplotype_evidence', [])
            
            if diplotype_evidence:
                builder.add_subsection_header("Diplotype-Specific Evidence")
                builder.add_spacer(0.08)
                
                from reportlab.platypus import Paragraph
                
                for idx, evidence in enumerate(diplotype_evidence[:5], 1):
                    drug = evidence.get('drug', 'Functional study')
                    category = evidence.get('category', '').upper()
                    
                    if category:
                        cat_color = self._get_category_color(category)
                        entry_header = f"<font color='{cat_color}'><b>[{category}]</b></font> {drug}"
                    else:
                        entry_header = drug
                    
                    builder.add_paragraph(entry_header, 'AlleleInfo')
                    
                    sentence = evidence.get('sentence', '').strip()
                    if sentence:
                        builder.add_paragraph(sentence, 'Small')
                    
                    metabolizer = evidence.get('metabolizer_types', '').strip()
                    if metabolizer:
                        builder.add_paragraph(f"<i>Metabolizer Type: {metabolizer}</i>", 'Small')
                    
                    pmid = evidence.get('pmid', '').strip()
                    if pmid and pmid not in ['nan', 'NaN', '']:
                        pmid_text = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>PMID: {pmid}</u></a>'
                        builder.add_paragraph(pmid_text, 'Small')
                    
                    if idx < len(diplotype_evidence[:5]):
                        builder.add_spacer(0.08)
                
                builder.add_spacer(0.15)
            
            # Drug interactions in table format - Grouped by therapeutic class
            if gene_entry['drug_recommendations'] or gene_entry['affected_drugs']:
                builder.add_subsection_header("Affected Drugs")
                builder.add_spacer(0.08)
                
                from reportlab.platypus import Paragraph
                
                # Deduplicate recommendations by drug name (keep first occurrence)
                seen_drugs = set()
                unique_recs = []
                for rec in gene_entry['drug_recommendations']:
                    drug = rec.get('drug', 'Unknown').title()
                    if drug not in seen_drugs:
                        seen_drugs.add(drug)
                        unique_recs.append(rec)
                
                # Sort by guideline level (highest priority first)
                unique_recs.sort(key=lambda x: self._get_guideline_level(x))
                
                # Group by therapeutic class
                drugs_by_class = {}
                for rec in unique_recs:
                    drug = rec.get('drug', 'Unknown').title()
                    therapeutic_class = self._classify_therapeutic_class(drug)
                    if therapeutic_class not in drugs_by_class:
                        drugs_by_class[therapeutic_class] = []
                    drugs_by_class[therapeutic_class].append(rec)
                
                # Show clinical recommendations table if available - by therapeutic class
                if unique_recs:
                    # Display drugs grouped by therapeutic class
                    for therapeutic_class in sorted(drugs_by_class.keys()):
                        class_recs = drugs_by_class[therapeutic_class]
                        
                        # Class header
                        builder.add_paragraph(f"<b><font color='#2E75B6'>{therapeutic_class}</font></b>", 'Normal')
                        builder.add_spacer(0.05)
                        
                        # Build comprehensive interaction table for this class
                        interaction_table = [
                            [
                                Paragraph('<b>Drug</b>', builder.styles['Normal']),
                                Paragraph('<b>Direction</b>', builder.styles['Normal']),
                                Paragraph('<b>Level</b>', builder.styles['Normal'])
                            ]
                        ]
                        
                        for rec in class_recs:
                            drug = rec.get('drug', 'Unknown').title()
                            direction = rec.get('direction', 'N/A')
                            guideline_level = self._get_guideline_level(rec)
                            level_label = ['A', 'B', 'C', 'D'][guideline_level - 1] if guideline_level <= 4 else 'D'
                            
                            # Build table row with Paragraph objects for proper formatting
                            interaction_table.append([
                                Paragraph(drug, builder.styles['Normal']),
                                Paragraph(direction or 'N/A', builder.styles['Normal']),
                                Paragraph(f"<b>{level_label}</b>", builder.styles['Normal'])
                            ])
                        
                        # Add table with proper column widths
                        col_widths = [2.0*inch, 3.5*inch, 0.6*inch]
                        builder.add_table(interaction_table, col_widths)
                        builder.add_spacer(0.12)
                
                # Add detailed explanations and clinical evidence below table
                if unique_recs:
                    builder.add_subsection_header("Clinical Evidence & Summary")
                    builder.add_spacer(0.08)
                    
                    for rec in unique_recs:
                        direction = rec.get('direction', '').strip()
                        reason = rec.get('reason', '').strip()
                        notes = rec.get('notes', '').strip()
                        pmid = rec.get('pmid', '').strip()
                        
                        # Only show explanation if there's actual evidence/direction
                        if direction and direction.lower() != 'unknown':
                            drug = rec.get('drug', 'Unknown').title()
                            
                            # Drug name as subsection
                            builder.add_paragraph(f"<b>{drug}</b>", 'AlleleInfo')
                            
                            # PMID reference with muted styling and clickable link
                            if pmid and pmid not in ['nan', 'NaN', '', 'None']:
                                try:
                                    pmid_num = int(pmid) if isinstance(pmid, str) else pmid
                                    pmid_link = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid_num}/" color="blue"><u>PMID: {pmid_num}</u></a>'
                                    builder.add_paragraph(pmid_link, 'Muted')
                                except (ValueError, TypeError):
                                    pass
                            
                            # Clinical summary narrative in smaller, denser text
                            if direction:
                                direction_lower = direction.lower()
                                
                                if 'decreased' in direction_lower:
                                    narrative = f"<b>Decreased metabolism.</b> This allele DECREASES how the body breaks down {drug}, meaning the patient needs a <b>LOWER dose</b>."
                                elif 'increased' in direction_lower:
                                    if 'toxicity' in direction_lower or 'adverse' in direction_lower or 'risk' in direction_lower:
                                        narrative = f"<b>Increased effect/toxicity.</b> This allele causes {drug} to have a STRONGER effect or increased adverse reactions. <b>AVOID or use with extreme caution</b>."
                                    else:
                                        narrative = f"<b>Increased metabolism.</b> This allele INCREASES how the body breaks down {drug}, meaning the patient may need a <b>HIGHER dose</b> for therapeutic effect."
                                else:
                                    narrative = f"<b>{direction.capitalize()}.</b> This allele affects the metabolism and efficacy of {drug} in the patient."
                                
                                builder.add_paragraph(narrative, 'Small')
                            
                            # Add detailed clinical notes if available (muted quote style)
                            if notes and notes not in ['nan', 'NaN', '']:
                                notes_formatted = notes.replace('\n', ' ')
                                builder.add_paragraph(f"“{notes_formatted}”", 'EvidenceNote')
                            
                            builder.add_spacer(0.10)
                
                # Show other affected drugs if they exist and haven't been covered
                rec_drugs_lower = {r.get('drug', '').lower() for r in unique_recs}
                other_drugs = [d for d in gene_entry['affected_drugs'] if d.lower() not in rec_drugs_lower]
                
                # Add gene-specific DDI pair heatmap if available
                if gene in gene_heatmap_paths:
                    heatmap_path = gene_heatmap_paths[gene]
                    if Path(heatmap_path).exists():
                        builder.add_spacer(0.1)
                        builder.add_paragraph(f"<b>Drug-Drug Interaction Heatmap</b>", 'SectionBreak')
                        builder.add_spacer(0.08)
                        try:
                            from reportlab.platypus import Image as RLImage
                            img = RLImage(heatmap_path, width=7.0*inch, height=5.5*inch)
                            builder.story.append(img)
                            builder.add_spacer(0.15)
                        except Exception as e:
                            print(f"Warning: Could not embed heatmap for {gene}: {e}")
                
                if other_drugs:
                    builder.add_paragraph(f"<b>Additional Drug Interactions ({len(other_drugs)}):</b>", 'SectionBreak')
                    builder.add_spacer(0.08)
                    
                    # Split comma-separated drugs and deduplicate
                    unique_other_drugs = set()
                    for drug_entry in other_drugs:
                        # Ensure drug_entry is a string and not empty
                        if not isinstance(drug_entry, str) or not drug_entry.strip():
                            continue
                        # Split by comma if present
                        drugs_list = [d.strip() for d in drug_entry.split(',')]
                        for drug in drugs_list:
                            # Only add non-empty, non-numeric drug names
                            if drug and not drug.isdigit() and drug.lower() not in rec_drugs_lower:
                                unique_other_drugs.add(drug)
                    
                    # Format as comma-separated paragraph
                    if unique_other_drugs:
                        drugs_text = ', '.join(sorted(unique_other_drugs))
                        builder.add_paragraph(drugs_text, 'Normal')
                    builder.add_spacer(0.15)
            
            # Show DrugBank drug-drug interactions for this gene's associated drugs
            drugbank_ddis = gene_entry.get('drugbank_ddis', [])
            if drugbank_ddis:
                # Limit to high and moderate severity for readability
                high_severity = [ddi for ddi in drugbank_ddis if ddi.get('severity') == 'HIGH']
                moderate_severity = [ddi for ddi in drugbank_ddis if ddi.get('severity') == 'MODERATE']
                
                # Show high severity first
                if high_severity:
                    builder.add_paragraph(f"<b>DrugBank High-Severity Drug Interactions ({len(high_severity)}):</b>", 'SectionBreak')
                    builder.add_spacer(0.08)
                    
                    ddi_table = [
                        [Paragraph('<b>Drug 1</b>', builder.styles['Normal']),
                         Paragraph('<b>Drug 2</b>', builder.styles['Normal']),
                         Paragraph('<b>Interaction</b>', builder.styles['Normal'])]
                    ]
                    
                    for ddi in high_severity[:20]:  # Limit to first 20
                        drug1 = ddi.get('drug_1', '')
                        drug2 = ddi.get('drug_2', '')
                        desc = ddi.get('description', '')[:200]  # Truncate long descriptions
                        
                        ddi_table.append([
                            Paragraph(drug1, builder.styles['Small']),
                            Paragraph(drug2, builder.styles['Small']),
                            Paragraph(desc, builder.styles['Small'])
                        ])
                    
                    col_widths = [1.3*inch, 1.3*inch, 3.5*inch]
                    builder.add_table(ddi_table, col_widths)
                    builder.add_spacer(0.12)
                
                # Show moderate severity summary
                if moderate_severity:
                    builder.add_paragraph(
                        f"<b>DrugBank Moderate-Severity Drug Interactions:</b> {len(moderate_severity)} interactions found "
                        f"(see full JSON report for complete details)", 
                        'Small'
                    )
                    builder.add_spacer(0.15)
            
            builder.add_spacer(0.2)
        
        builder.build()
    
    def _setup_interaction_styles(self, builder: PDFReportBuilder) -> None:
        """Setup custom styles for drug interaction formatting"""
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
        
        styles = builder.styles
        
        # Centered gene header style
        if 'CenteredGeneHeader' not in styles:
            styles.add(ParagraphStyle(
                name='CenteredGeneHeader',
                parent=styles['Heading2'],
                alignment=TA_CENTER,
                spaceAfter=6,
                spaceBefore=12,
                fontSize=14,
                textColor='#003366'
            ))
        
        # Interaction-specific paragraph style
        if 'InteractionNarrative' not in styles:
            styles.add(ParagraphStyle(
                name='InteractionNarrative',
                parent=styles['Normal'],
                fontSize=9,
                spaceAfter=8,
                spaceBefore=2,
                alignment=TA_LEFT,
                leading=12
            ))
    
    def _add_gene_section_to_pdf(self, builder: PDFReportBuilder, gene: str, gene_data: Dict[str, Any]):
        """Add a gene section to the PDF report with proper spacing and page breaks"""
        
        # Gene header
        builder.add_spacer(0.12)
        builder.add_paragraph(f"<u><b>{gene}</b></u>", 'GeneHeader')
        
        # Variants section with monospace font and consequence type
        builder.add_paragraph(f"<b>Variants Detected ({gene_data['total_variants']}):</b>", 'SectionBreak')
        builder.add_spacer(0.08)
        
        for var in gene_data.get('variants', [])[:3]:
            # Format variant in monospace font
            var_notation = f"<font face='Courier'>{var['position']} {var['ref']}→{var['alt']}</font>"
            consequence = var.get('consequence', 'unknown')
            rsid = var.get('rsid', '-')
            
            var_text = f"• {var_notation} | <b>{consequence}</b> | rsID: {rsid}"
            builder.add_paragraph(var_text, 'Normal')
        
        builder.add_spacer(0.15)
        
        # Identified star alleles section
        # Only show the two alleles that make up the diplotype (from the diplotype call)
        diplotype_alleles = []
        if gene_data.get('diplotype'):
            # Extract the two alleles from diplotype format like "allele1/allele2"
            diplotype_alleles = gene_data['diplotype'].split('/')
        
        if diplotype_alleles:
            builder.add_paragraph("<b>IDENTIFIED STAR ALLELES & FUNCTIONALITY:</b>", 'SectionBreak')
            builder.add_spacer(0.1)
            
            for allele in sorted(diplotype_alleles):
                func_info = gene_data.get('allele_functionality', {}).get(allele, {})
                
                # Allele header with function status
                function_status = func_info.get('clinical_function', 'Unknown Function')
                evidence_level = func_info.get('strength_of_evidence', '')
                
                # Build compact header
                if evidence_level and evidence_level.lower() not in ['no evidence', 'none', 'unknown']:
                    allele_header = f"<b>{allele}: {function_status}</b> (Evidence: {evidence_level})"
                else:
                    allele_header = f"<b>{allele}: {function_status}</b>"
                
                builder.add_paragraph(allele_header, 'GeneHeader')
                
                # Only show summary if there's meaningful information (not just "unknown function" boilerplate)
                summary = func_info.get('summary', '').strip()
                if summary:
                    # Skip generic "unknown function" summaries that don't add value
                    skip_phrases = [
                        'assigned unknown function due to no literature',
                        'consensus among experts was unknown function',
                        'no available evidence'
                    ]
                    summary_lower = summary.lower()
                    is_generic = any(phrase in summary_lower for phrase in skip_phrases)
                    
                    # Only show summary if it contains actual information
                    if not is_generic or len(summary) > 150:  # Show if detailed enough
                        builder.add_spacer(0.05)
                        builder.add_paragraph(f"<b>Summary:</b> {summary}", 'Normal')
                
                builder.add_spacer(0.15)
        
        builder.add_spacer(0.15)
        
        # Diplotype prediction section
        if gene_data.get('diplotype'):
            builder.add_paragraph("<b>DIPLOTYPE PREDICTION:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            
            # Build diplotype prediction table
            table_data = [
                ['Property', 'Value']
            ]
            
            table_data.append(['Genotype', gene_data['diplotype']])
            
            if gene_data.get('coded_summary') and gene_data.get('coded_summary') != 'Unknown':
                table_data.append(['Phenotype', gene_data['coded_summary']])
            
            if gene_data.get('ehr_priority') and gene_data.get('ehr_priority') != 'Unknown':
                table_data.append(['EHR Priority', gene_data.get('ehr_priority')])
            
            if gene_data.get('metabolizer_status'):
                table_data.append(['Metabolizer Status', gene_data['metabolizer_status']])
            
            # Add table with proper styling
            col_widths = [1.5*inch, 4.5*inch]
            builder.add_table(table_data, col_widths)
            
            builder.add_spacer(0.15)
            
            # Add diplotype-specific evidence from PharmGKB
            diplotype_evidence = gene_data.get('diplotype_evidence', [])
            
            if diplotype_evidence:
                builder.add_paragraph("<b>Diplotype-Specific Evidence from PharmGKB:</b>", 'SectionBreak')
                builder.add_spacer(0.08)
                
                for idx, evidence in enumerate(diplotype_evidence[:5], 1):
                    drug = evidence.get('drug', 'Unknown drug')
                    category = evidence.get('category', '').upper()
                    
                    if category:
                        cat_color = self._get_category_color(category)
                        entry_header = f"<font color='{cat_color}'><b>[{category}]</b></font> <b>{drug}</b>"
                    else:
                        entry_header = f"<b>{drug}</b>"
                    
                    builder.add_paragraph(entry_header, 'Normal')
                    
                    sentence = evidence.get('sentence', '').strip()
                    if sentence:
                        builder.add_paragraph(sentence, 'Small')
                    
                    metabolizer = evidence.get('metabolizer_types', '').strip()
                    if metabolizer:
                        builder.add_paragraph(f"<i>Metabolizer Type: {metabolizer}</i>", 'Small')
                    
                    pmid = evidence.get('pmid', '').strip()
                    if pmid and pmid not in ['nan', 'NaN', '']:
                        pmid_text = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>PMID: {pmid}</u></a>'
                        builder.add_paragraph(pmid_text, 'Small')
                    
                    if idx < len(diplotype_evidence[:5]):
                        builder.add_spacer(0.08)
            else:
                builder.add_paragraph(
                    f"<i>No diplotype-specific evidence found for {gene} {gene_data['diplotype']} in PharmGKB databases.</i>",
                    'Normal'
                )
                builder.add_spacer(0.08)
        
        builder.add_spacer(0.15)
        
        # Affected drugs summary
        if gene_data.get('affected_drugs'):
            affected_count = len(gene_data['affected_drugs'])
            spec_recs = gene_data.get('pharmgkb_recommendations', [])
            
            if not spec_recs:
                # No specific interactions found for this diplotype
                builder.add_paragraph(
                    f"<b>Drug Interactions:</b> No diplotype specific drug is found on pharmGKB. Total {affected_count} drug(s) documented.",
                    'SectionBreak'
                )
                builder.add_spacer(0.1)
        
        # Clinical Evidence & Supporting Literature - Grouped by Strength
        if gene_data.get('clinical_evidence'):
            builder.add_paragraph("<b>CLINICAL EVIDENCE & SUPPORTING LITERATURE:</b>", 'SectionBreak')
            builder.add_spacer(0.08)
            
            # Group evidence by strength
            evidence_by_strength = {'Level A': [], 'Level B': [], 'Level C': []}
            for evidence in gene_data['clinical_evidence']:
                strength = self._classify_evidence_strength(evidence)
                evidence_by_strength[strength].append(evidence)
            
            # Display evidence grouped by strength
            global_idx = 1
            for strength_level in ['Level A', 'Level B', 'Level C']:
                evidence_list = evidence_by_strength[strength_level]
                if not evidence_list:
                    continue
                
                # Strength level header
                strength_colors = {
                    'Level A': '#C0392B',
                    'Level B': '#F39C12', 
                    'Level C': '#27AE60'
                }
                color = strength_colors.get(strength_level, '#7F8C8D')
                builder.add_paragraph(f"<font color='{color}'><b>{strength_level} Evidence ({len(evidence_list)} items):</b></font>", 'Normal')
                builder.add_spacer(0.08)
                
                evidence_rows = [
                    [
                        Paragraph('<b>#</b>', builder.styles['Small']),
                        Paragraph('<b>Context</b>', builder.styles['Small']),
                        Paragraph('<b>Summary</b>', builder.styles['Small']),
                        Paragraph('<b>PMID</b>', builder.styles['Small'])
                    ]
                ]

                for evidence in evidence_list:
                    evidence_type = evidence.get('type', 'Unknown').upper()
                    type_color = self._get_evidence_type_color(evidence_type)
                    allele_raw = evidence.get('allele', '')
                    # Filter allele to only show rsIDs (remove COSV, CM, etc.)
                    allele_parts = [p.strip() for p in allele_raw.replace('&', ';').split(';')]
                    allele = ';'.join([p for p in allele_parts if p.startswith('rs')]) if any(p.startswith('rs') for p in allele_parts) else allele_raw
                    category = evidence.get('category', '')

                    # Build context with category badge and drug only
                    context_parts = []
                    
                    # Add category badge if available
                    if category:
                        cat_color = self._get_category_color(category)
                        context_parts.append(f"<font color='{cat_color}'><b>[{category.upper()}]</b></font>")
                    
                    # Add drug in bold - this is the only descriptive text we show
                    if evidence.get('drug'):
                        context_parts.append(f"<b>{evidence['drug']}</b>")
                    
                    # Format context with line breaks for better readability
                    if context_parts:
                        context_str = "<br/>".join(context_parts)
                    else:
                        context_str = "<i>Functional study</i>"

                    # Use Sentence field if available (cleaner summary), otherwise use notes
                    sentence = (evidence.get('sentence') or '').strip()
                    notes = (evidence.get('notes') or '').strip()
                    
                    if sentence and sentence not in ['nan', 'NaN', '']:
                        summary_txt = sentence.replace('\n', ' ').replace('\r', ' ')
                    elif notes and notes not in ['nan', 'NaN']:
                        summary_txt = notes.replace('\n', ' ').replace('\r', ' ')
                    else:
                        summary_txt = '<i>No summary provided</i>'

                    # Make PMID clickable
                    pmid = evidence.get('pmid', '').strip()
                    if pmid and pmid not in ['nan', 'NaN', '']:
                        pmid_display = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" color="blue"><u>{pmid}</u></a>'
                    else:
                        pmid_display = '—'

                    evidence_rows.append([
                        Paragraph(f"<font color='{type_color}'><b>{global_idx}</b></font>", builder.styles['Small']),
                        Paragraph(f"<b>{allele}</b><br/><font color='{type_color}'><b>{evidence_type}</b></font><br/>{context_str}", builder.styles['Small']),
                        Paragraph(summary_txt, builder.styles['Small']),
                        Paragraph(pmid_display, builder.styles['Small'])
                    ])
                    global_idx += 1

                # Add evidence table for this strength level
                col_widths = [0.4*inch, 2.1*inch, 2.5*inch, 0.9*inch]
                builder.add_table(evidence_rows, col_widths)
                builder.add_spacer(0.15)

            builder.add_spacer(0.1)
        
        # Clinical notes
        
        if gene_data.get('clinical_notes'):
            builder.add_spacer(0.1)
            for note in gene_data['clinical_notes']:
                builder.add_paragraph(f"⚠ <b>Note:</b> {note}", 'Normal')
        
        builder.add_spacer(0.3)
    
    def _get_evidence_type_color(self, evidence_type: str) -> str:
        """Get color for evidence type"""
        colors = {
            'PHENOTYPE': '#E74C3C',
            'DRUG': '#3498DB',
            'FUNCTIONAL': '#2ECC71'
        }
        return colors.get(evidence_type, '#95A5A6')
    
    def _get_category_color(self, category: str) -> str:
        """Get color for phenotype category badges"""
        category_upper = category.upper()
        colors = {
            'DOSAGE': '#2980B9',
            'METABOLISM/PK': '#27AE60',
            'TOXICITY': '#C0392B',
            'EFFICACY': '#8E44AD',
            'OTHER': '#7F8C8D'
        }
        return colors.get(category_upper, '#34495E')
    
    def _get_priority_color(self, priority: str) -> str:
        """Get color for EHR priority level"""
        priority_upper = priority.upper() if priority else 'NONE'
        if priority_upper in ['HIGH', 'REQUIRED']:
            return '#C0392B'  # Red
        elif priority_upper in ['MEDIUM', 'MODERATE']:
            return '#F39C12'  # Yellow/Orange
        elif priority_upper in ['LOW', 'NORMAL', 'ROUTINE']:
            return '#27AE60'  # Green
        else:
            return '#7F8C8D'  # Gray
    
    def _classify_therapeutic_class(self, drug_name: str) -> str:
        """Classify drugs into therapeutic categories"""
        drug_lower = drug_name.lower()
        
        # Cardiovascular
        if any(x in drug_lower for x in ['warfarin', 'clopidogrel', 'atorvastatin', 'simvastatin', 'pravastatin']):
            return 'Cardiovascular'
        
        # Psychiatric/CNS
        if any(x in drug_lower for x in ['citalopram', 'escitalopram', 'sertraline', 'fluoxetine', 'paroxetine', 
                                          'venlafaxine', 'amitriptyline', 'nortriptyline', 'desipramine', 'imipramine',
                                          'clomipramine', 'doxepin', 'haloperidol', 'risperidone', 'aripiprazole',
                                          'olanzapine', 'clozapine', 'codeine', 'tramadol', 'oxycodone', 'hydrocodone']):
            return 'Psychiatric/CNS'
        
        # Oncology
        if any(x in drug_lower for x in ['fluorouracil', '5-fu', 'capecitabine', 'azathioprine', 'mercaptopurine',
                                          'thioguanine', 'irinotecan', 'tamoxifen', 'cisplatin']):
            return 'Oncology'
        
        # Anticoagulation
        if any(x in drug_lower for x in ['warfarin', 'acenocoumarol', 'phenprocoumon']):
            return 'Anticoagulation'
        
        # Pain Management
        if any(x in drug_lower for x in ['codeine', 'tramadol', 'oxycodone', 'hydrocodone', 'morphine', 'fentanyl']):
            return 'Pain Management'
        
        # Infectious Disease
        if any(x in drug_lower for x in ['abacavir', 'efavirenz', 'atazanavir', 'isoniazid', 'rifampin']):
            return 'Infectious Disease'
        
        # Immunosuppression
        if any(x in drug_lower for x in ['azathioprine', 'mercaptopurine', 'tacrolimus', 'cyclosporine']):
            return 'Immunosuppression'
        
        return 'Other'
    
    def _get_guideline_level(self, recommendation: Dict[str, Any]) -> int:
        """Get CPIC/PharmGKB guideline level for sorting (lower number = higher priority)"""
        # Check for explicit level indicators
        direction = str(recommendation.get('direction', '')).upper()
        category = str(recommendation.get('category', '')).upper()
        
        # Level A (highest priority): Strong recommendations, dosing changes, avoid
        if any(x in direction for x in ['AVOID', 'CONTRAINDICATED', 'ALTERNATIVE REQUIRED']):
            return 1
        if 'DOSAGE' in category or 'DOSE' in direction:
            return 2
        if 'TOXICITY' in category or 'ADVERSE' in direction:
            return 2
        
        # Level B: Moderate recommendations
        if 'MONITOR' in direction or 'CAUTION' in direction:
            return 3
        if 'EFFICACY' in category:
            return 3
        
        # Level C/D: Lower priority
        return 4
    
    def _classify_evidence_strength(self, evidence: Dict[str, Any]) -> str:
        """Classify evidence by strength (Level A, B, C)"""
        # Check for explicit evidence level
        evidence_level = evidence.get('level_of_evidence', '').upper()
        if evidence_level in ['1A', '1B', 'A', 'HIGH']:
            return 'Level A'
        elif evidence_level in ['2A', '2B', 'B', 'MODERATE']:
            return 'Level B'
        elif evidence_level in ['3', '4', 'C', 'LOW', 'D']:
            return 'Level C'
        
        # Infer from evidence type and source
        evidence_type = evidence.get('type', '').upper()
        category = evidence.get('category', '').upper()
        
        # High quality: Clinical studies with dosage/toxicity implications
        if category in ['DOSAGE', 'TOXICITY'] and 'DRUG' in evidence_type:
            return 'Level A'
        
        # Moderate quality: Efficacy studies, metabolism studies
        if category in ['EFFICACY', 'METABOLISM/PK']:
            return 'Level B'
        
        # Lower quality: Functional studies, phenotype associations
        if 'FUNCTIONAL' in evidence_type or 'PHENOTYPE' in evidence_type:
            return 'Level C'
        
        return 'Level C'
    
    def _generate_executive_summary_page(self, builder: PDFReportBuilder, analysis_results: Dict[str, Any]):
        """Generate executive summary page with key findings"""
        from reportlab.platypus import Paragraph, Spacer
        from reportlab.lib.units import inch
        from reportlab.lib import colors
        
        # Title
        builder.add_paragraph("<b>EXECUTIVE SUMMARY</b>", 'Heading1')
        builder.add_spacer(0.2)
        
        # Patient metadata
        metadata_text = f"""
        <b>Sample ID:</b> {analysis_results['sample_id']}<br/>
        <b>Analysis Date:</b> {analysis_results['analysis_date']}<br/>
        <b>Analysis ID:</b> {analysis_results['analysis_id']}
        """
        builder.add_paragraph(metadata_text, 'Normal')
        builder.add_spacer(0.2)
        
        # Introduction
        intro_text = f"""This pharmacogenomics report analyzes <b>{len(analysis_results['genes_analyzed'])} genes</b> to identify genetic variants that may affect drug metabolism and response."""
        builder.add_paragraph(intro_text, 'Normal')
        builder.add_spacer(0.25)
        
        # Classify genes by priority
        high_priority_genes = []
        medium_priority_genes = []
        normal_priority_genes = []
        
        for gene in sorted(analysis_results['variants_by_gene'].keys()):
            gene_data = analysis_results['variants_by_gene'][gene]
            priority = gene_data.get('ehr_priority', 'none')
            
            if priority and isinstance(priority, str):
                priority_upper = priority.upper()
                if priority_upper in ['HIGH', 'REQUIRED']:
                    high_priority_genes.append(gene)
                elif priority_upper in ['MEDIUM', 'MODERATE']:
                    medium_priority_genes.append(gene)
                elif priority_upper in ['LOW', 'NORMAL', 'ROUTINE']:
                    normal_priority_genes.append(gene)
                else:
                    normal_priority_genes.append(gene)
            else:
                normal_priority_genes.append(gene)
        
        # Key Findings header
        builder.add_paragraph("<b>Key Findings:</b>", 'Heading3')
        builder.add_spacer(0.12)
        
        # Build findings list with color-coded priorities
        findings = []
        if high_priority_genes:
            count = len(high_priority_genes)
            findings.append(f"<font color='#C0392B'><b>• {count} gene{'s' if count != 1 else ''} with HIGH-PRIORITY variants requiring clinical attention</b></font>")
        
        if medium_priority_genes:
            count = len(medium_priority_genes)
            findings.append(f"<font color='#F39C12'><b>• {count} gene{'s' if count != 1 else ''} with MEDIUM-PRIORITY variants</b></font>")
        
        if normal_priority_genes:
            count = len(normal_priority_genes)
            findings.append(f"<font color='#27AE60'>• {count} gene{'s' if count != 1 else ''} with normal or routine variants</font>")
        
        for finding in findings:
            builder.add_paragraph(finding, 'Normal')
            builder.add_spacer(0.08)
        
        builder.add_spacer(0.2)
        
        # High Priority Genes detail
        if high_priority_genes:
            builder.add_paragraph("<b><font color='#C0392B'>High Priority Genes:</font></b>", 'Heading3')
            builder.add_spacer(0.1)
            genes_text = ", ".join(sorted(high_priority_genes))
            builder.add_paragraph(f"<font color='#C0392B'><b>{genes_text}</b></font>", 'Normal')
            builder.add_spacer(0.15)
        
        if medium_priority_genes:
            builder.add_paragraph("<b><font color='#F39C12'>Medium Priority Genes:</font></b>", 'Heading3')
            builder.add_spacer(0.1)
            genes_text = ", ".join(sorted(medium_priority_genes))
            builder.add_paragraph(f"<font color='#F39C12'><b>{genes_text}</b></font>", 'Normal')
            builder.add_spacer(0.15)
        
        # Add disclaimer
        builder.add_spacer(0.3)
        disclaimer = """<i><b>Note:</b> This report provides genetic information for clinical decision support. All treatment decisions should be made by qualified healthcare professionals considering the complete clinical context.</i>"""
        builder.add_paragraph(disclaimer, 'Small')
        
        # Page break after executive summary
        builder.add_page_break()
                    

def main():
    """Main analysis function"""
    
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Pharmacogenomics Analysis Pipeline')
    parser.add_argument('--input', type=str, help='Path to input annotation TSV file')
    parser.add_argument('--source', type=str, help='Path to source data directory (contains source-allele-definition, pharmgkb, etc.)')
    parser.add_argument('--output', type=str, help='Path to output directory for reports')
    parser.add_argument('--drug-db', type=str, help='Path to drug interaction database (txt or DrugBank XML)')
    parser.add_argument('--meds', type=str, help='Comma-separated list of patient medications (e.g. "Warfarin,Simvastatin")')
    
    args = parser.parse_args()
    
    # Resolve to absolute path of script location
    script_dir = Path(__file__).resolve().parent
    
    # Determine input file
    if args.input:
        sample_file = Path(args.input)
        if not sample_file.exists():
            print(f"Error: Input file not found: {sample_file}")
            return
    else:
        # Default to sample_annotation.tsv in script directory or parent
        sample_file = script_dir / "sample_annotation.tsv"
        if not sample_file.exists():
            parent_sample = script_dir.parent / "sample_annotation.tsv"
            if parent_sample.exists():
                sample_file = parent_sample
            else:
                print(f"Error: Sample file not found in {sample_file} or {parent_sample}")
                print(f"Script directory: {script_dir}")
                return
    
    # Determine source directory
    if args.source:
        source_dir = Path(args.source)
        if not source_dir.exists():
            print(f"Error: Source directory not found: {source_dir}")
            return
    else:
        # Default to parent of script directory (pharmacogenomics/)
        source_dir = sample_file.parent if sample_file.parent.name == "pharmacogenomics" else script_dir.parent
    
    # Determine output directory
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        # Default to source directory
        output_dir = source_dir
    
    try:
        # Initialize analyzer
        analyzer = ImprovedPharmacogenomicsAnalyzer(
            sample_annotation_file=str(sample_file),
            source_dir=str(source_dir)
        )
        
        # Override output directory
        analyzer.source_dir = output_dir
        
        # Run analysis
        print("Starting improved pharmacogenomics analysis...\n")
        results = analyzer.analyze_sample()
        print(f"\n✓ Variant analysis complete. Analyzing {len(results.get('genes_analyzed', []))} genes")
        
        # Optionally integrate drug interaction database (text dump or DrugBank XML)
        if args.drug_db:
            print(f"\n📋 Drug interaction database specified: {args.drug_db}")
            db_path = Path(args.drug_db)
            if not db_path.exists():
                print(f"Warning: Drug DB file not found at {db_path}. Skipping interactions.")
                patient_meds = []  # Initialize for Phase 2
            else:
                print(f"✓ Found DrugBank XML at: {db_path}")
                # Determine patient meds: CLI-provided meds take precedence; otherwise derive from PharmGKB annotations
                if args.meds and args.meds.strip().upper() != 'NONE':
                    patient_meds = [m.strip() for m in args.meds.split(',') if m.strip()]
                else:
                    print("   Deriving patient medications from analysis results...")
                    patient_meds_set = set()
                    # Derive medications from the analysis results                    
                    for gene, gene_data in results.get('variants_by_gene', {}).items():
                        # Add from affected_drugs field
                        for drug in gene_data.get('affected_drugs', []):
                            if drug and drug.strip():
                                patient_meds_set.add(drug.strip())
                        
                        # Add from pharmgkb_recommendations
                        for rec in gene_data.get('pharmgkb_recommendations', []):
                            drug = rec.get('drug', '')
                            if drug and drug.strip():
                                patient_meds_set.add(drug.strip())
                    
                    # Also try PharmGKB cache if available
                    if analyzer.pharmgkb_cache:
                        for gene in results.get('variants_by_gene', {}).keys():
                            try:
                                drugs = analyzer.pharmgkb_cache.get_affected_drugs(gene)
                                for d in drugs:
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
                    # For XML, pass the path string so integrator can stream-parse; otherwise read text
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
            patient_meds = []  # Initialize for Phase 2
        
        # Integrate Phase 2 enhancements
        print("\n🔬 Integrating Phase 2 analysis (DPWG, FDA warnings, PD genes, dosing guidance)...")
        try:
            results = analyzer.integrate_phase2_analysis(results, patient_meds)
            print("   ✓ Phase 2 integration complete")
            
            # Print Phase 2 summary
            phase2_summary = results.get('phase2_analysis', {}).get('phase2_summary', {})
            if phase2_summary:
                print(f"\n📊 Phase 2 Summary:")
                print(f"   Non-star genes analyzed: {phase2_summary.get('total_non_star_genes_analyzed', 0)}")
                print(f"   DPWG recommendations: {phase2_summary.get('dpwg_recommendations_found', 0)}")
                print(f"   FDA warnings: {phase2_summary.get('fda_warnings_triggered', 0)}")
                print(f"   Dosing adjustments: {phase2_summary.get('dosing_adjustments_needed', 0)}")
        except Exception as e:
            print(f"   ⚠️  Phase 2 integration error: {e}")
            import traceback
            traceback.print_exc()
        
        # Print priority summary
        print("\n📊 Generating priority summary...")
        analyzer._print_priority_summary(results)
        
        # Generate reports
        analyzer.generate_reports(results)
        
        print("\n✓ Analysis complete!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
