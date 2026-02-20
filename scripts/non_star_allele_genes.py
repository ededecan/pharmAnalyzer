#!/usr/bin/env python3
"""
Non-Star Allele Gene Handler
Handles genes that don't use star allele nomenclature (VKORC1, CYP4F2, IFNL3, etc.)
- Reports raw variant information
- Includes PharmGKB annotations
- Uses IMPACT for functional predictions
- Extracts drug-specific recommendations
"""

from typing import Dict, List, Optional, Any, Set, Tuple
from dataclasses import dataclass
from collections import defaultdict
import re


@dataclass
class VariantAnnotation:
    """Raw variant with PharmGKB annotations"""
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    rsid: Optional[str]
    impact: str  # HIGH, MODERATE, LOW
    consequence: str
    allele_frequency: Optional[float]
    genotype: Optional[str]  # e.g., "C/C", "C/T", "T/T"
    

@dataclass
class PharmGKBVariantAnnotation:
    """PharmGKB annotation for a variant"""
    variant_id: str
    drug: str
    association: str  # "Associated", "Not associated"
    direction: str  # "increased", "decreased", "no change"
    effect_type: str  # "dose", "response", "metabolism", "efficacy"
    evidence_level: str  # "A", "B", "C"
    phenotype_category: str
    pmid: str
    sentence: str
    

@dataclass  
class NonStarAlleleFinding:
    """Finding for non-star-allele gene"""
    gene: str
    variants: List[VariantAnnotation]
    phenotype_predictions: Dict[str, Any]
    drug_recommendations: Dict[str, List[Dict[str, Any]]]
    warfarin_dosing_guidance: Optional[str]  # For VKORC1
    hepatitis_response: Optional[str]  # For IFNL3
    clinical_significance: str  # HIGH, MODERATE, LOW
    

class NonStarAlleleGeneHandler:
    """
    Handles analysis of genes without star allele nomenclature
    Primary targets: VKORC1, CYP4F2, IFNL3
    """
    
    # Genes that don't use star allele nomenclature
    NON_STAR_ALLELE_GENES = {
        'VKORC1', 'CYP4F2', 'IFNL3', 'SLCO1B1', 
        'F2', 'F5', 'MTHFR', 'PAH', 'NUDT15'
    }
    
    # Impact severity mapping
    IMPACT_PRIORITY = {
        'HIGH': 0,
        'MODERATE': 1,
        'LOW': 2,
        'MODIFIER': 3
    }
    
    # VKORC1 variant-specific warfarin dosing information
    VKORC1_VARIANTS = {
        'rs9923231': {
            'genotype_effects': {
                'GG': {'description': '-1639G/G', 'warfarin_sensitivity': 'high', 'dose_prediction': 'lower'},
                'GA': {'description': '-1639G/A', 'warfarin_sensitivity': 'intermediate', 'dose_prediction': 'intermediate'},
                'AA': {'description': '-1639A/A', 'warfarin_sensitivity': 'lower', 'dose_prediction': 'higher'}
            },
            'gene_description': 'Vitamin K epoxide reductase complex subunit 1',
            'clinical_relevance': 'Major warfarin dosing determinant',
            'reference_url': 'https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/'
        }
    }
    
    # IFNL3 (IL28B) hepatitis C response
    IFNL3_VARIANTS = {
        'rs12979860': {
            'genotype_effects': {
                'CC': {'description': 'IL28B CC', 'response_rate': '70-80%', 'recommendation': 'Standard therapy'},
                'CT': {'description': 'IL28B CT', 'response_rate': '40-50%', 'recommendation': 'Consider extended therapy'},
                'TT': {'description': 'IL28B TT', 'response_rate': '25-30%', 'recommendation': 'Extended therapy or DAA'}
            },
            'gene_description': 'Interferon-lambda 3',
            'clinical_relevance': 'Hepatitis C spontaneous clearance and treatment response',
            'reference_url': 'https://www.ncbi.nlm.nih.gov/grc/human/genes/282617'
        },
        'rs8099917': {
            'genotype_effects': {
                'TT': {'description': 'IL28B TT', 'response_rate': '70-80%', 'recommendation': 'Standard therapy'},
                'TG': {'description': 'IL28B TG', 'response_rate': '40-50%', 'recommendation': 'Consider extended therapy'},
                'GG': {'description': 'IL28B GG', 'response_rate': '25-30%', 'recommendation': 'Extended therapy or DAA'}
            },
            'gene_description': 'Interferon-lambda 3',
            'clinical_relevance': 'Hepatitis C spontaneous clearance and treatment response',
            'reference_url': 'https://www.ncbi.nlm.nih.gov/grc/human/genes/282617'
        }
    }
    
    # CYP4F2 warfarin sensitivity (secondary to VKORC1/CYP2C9)
    CYP4F2_VARIANTS = {
        'rs2108622': {
            'genotype_effects': {
                'CC': {'description': 'CYP4F2 *2/*2 (normal)', 'warfarin_effect': 'baseline'},
                'CT': {'description': 'CYP4F2 *1/*2', 'warfarin_effect': 'slight increase'},
                'TT': {'description': 'CYP4F2 *1/*1', 'warfarin_effect': 'increased dose requirement'}
            },
            'gene_description': 'Cytochrome P450 family 4 subfamily F member 2',
            'clinical_relevance': 'Secondary warfarin dosing factor (~5% explained variance)',
            'reference_url': 'https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/'
        }
    }
    
    # Drug recommendations for non-star genes
    DRUG_RECOMMENDATIONS = {
        'VKORC1': {
            'warfarin': {
                'high_sensitivity': 'Use lower starting dose (2-3mg daily), more frequent INR monitoring',
                'intermediate_sensitivity': 'Use standard dosing with pharmacogenomic prediction',
                'lower_sensitivity': 'May require higher doses, standard monitoring'
            },
            'acenocoumarol': 'Consider pharmacogenomic dosing prediction',
            'phenprocoumon': 'Consider pharmacogenomic dosing prediction'
        },
        'IFNL3': {
            'interferon_alpha': 'CC genotype: higher cure rates with IFN/RBV',
            'ribavirin': 'CC genotype: improved sustained virologic response',
            'sofosbuvir': 'Effective for all genotypes; DAA preferred for TT carriers',
            'ledipasvir_sofosbuvir': 'Recommended for all IFNL3 genotypes'
        },
        'CYP4F2': {
            'warfarin': 'Secondary dosing factor; use with VKORC1 and CYP2C9 results'
        }
    }
    
    def __init__(self):
        """Initialize the handler"""
        self.pharmgkb_annotations: Dict[str, List[PharmGKBVariantAnnotation]] = {}
        
    def is_non_star_allele_gene(self, gene: str) -> bool:
        """Check if gene uses non-star-allele nomenclature"""
        return gene in self.NON_STAR_ALLELE_GENES
    
    def parse_variant_info(self, variant_record: Dict[str, Any]) -> VariantAnnotation:
        """
        Parse variant information from VCF/annotation record
        
        Args:
            variant_record: Dictionary with variant info
                - chromosome: str
                - position: int
                - ref_allele: str
                - alt_allele: str
                - rsid: str (optional)
                - impact: str (HIGH/MODERATE/LOW)
                - consequence: str
                - allele_frequency: float (optional)
                - genotype: str (e.g., "C/T")
                
        Returns:
            VariantAnnotation object
        """
        return VariantAnnotation(
            chromosome=variant_record.get('chromosome'),
            position=variant_record.get('position'),
            ref_allele=variant_record.get('ref_allele'),
            alt_allele=variant_record.get('alt_allele'),
            rsid=variant_record.get('rsid'),
            impact=variant_record.get('impact', 'MODIFIER'),
            consequence=variant_record.get('consequence', 'Unknown'),
            allele_frequency=variant_record.get('allele_frequency'),
            genotype=variant_record.get('genotype')
        )
    
    def get_vkorc1_warfarin_guidance(self, variants: List[VariantAnnotation]) -> str:
        """
        Generate warfarin dosing guidance for VKORC1 rs9923231 variant
        
        Args:
            variants: List of VariantAnnotation objects for VKORC1
            
        Returns:
            Human-readable dosing guidance
        """
        guidance_lines = ["WARFARIN DOSING GUIDANCE (VKORC1):\n"]
        
        for var in variants:
            if var.rsid == 'rs9923231' or '-1639' in var.consequence:
                if var.genotype:
                    variant_key = 'rs9923231'
                    if variant_key in self.VKORC1_VARIANTS:
                        info = self.VKORC1_VARIANTS[variant_key]
                        if var.genotype in info['genotype_effects']:
                            effect = info['genotype_effects'][var.genotype]
                            guidance_lines.append(f"  Genotype {var.genotype} ({effect['description']}):")
                            guidance_lines.append(f"    - Warfarin Sensitivity: {effect['warfarin_sensitivity'].upper()}")
                            guidance_lines.append(f"    - Dose Prediction: {effect['dose_prediction'].upper()} doses")
                            guidance_lines.append(f"    - Expected INR: Monitor closely (target 2-3)")
        
        if len(guidance_lines) == 1:
            guidance_lines.append("  Unable to determine - use standard warfarin dosing")
        
        return "\n".join(guidance_lines)
    
    def get_ifnl3_hcv_response(self, variants: List[VariantAnnotation]) -> str:
        """
        Generate HCV treatment response prediction for IFNL3
        
        Args:
            variants: List of VariantAnnotation objects for IFNL3
            
        Returns:
            Human-readable HCV response guidance
        """
        guidance_lines = ["HEPATITIS C TREATMENT RESPONSE (IFNL3):\n"]
        
        for var in variants:
            if var.rsid in self.IFNL3_VARIANTS:
                info = self.IFNL3_VARIANTS[var.rsid]
                if var.genotype and var.genotype in info['genotype_effects']:
                    effect = info['genotype_effects'][var.genotype]
                    guidance_lines.append(f"  Variant {var.rsid} - Genotype {var.genotype}:")
                    guidance_lines.append(f"    - Phenotype: {effect['description']}")
                    guidance_lines.append(f"    - HCV Response Rate: {effect['response_rate']}")
                    guidance_lines.append(f"    - Recommendation: {effect['recommendation']}")
        
        if len(guidance_lines) == 1:
            guidance_lines.append("  No major IFNL3/IL28B variants detected")
        
        return "\n".join(guidance_lines)
    
    def analyze_gene_variants(self, 
                             gene: str, 
                             variants: List[VariantAnnotation],
                             pharmgkb_annotations: Dict[str, List[str]] = None) -> NonStarAlleleFinding:
        """
        Comprehensive analysis of non-star-allele gene
        
        Args:
            gene: Gene name (VKORC1, CYP4F2, IFNL3, etc.)
            variants: List of detected variants
            pharmgkb_annotations: PharmGKB drug annotations (optional)
            
        Returns:
            NonStarAlleleFinding with complete analysis
        """
        # Sort variants by impact
        sorted_variants = sorted(
            variants,
            key=lambda v: self.IMPACT_PRIORITY.get(v.impact, 99)
        )
        
        # Determine clinical significance
        clinical_significance = self._determine_clinical_significance(gene, sorted_variants)
        
        # Get gene-specific recommendations
        drug_recommendations = self._extract_drug_recommendations(gene, sorted_variants)
        
        # Get specialized guidance
        warfarin_guidance = None
        hcv_guidance = None
        
        if gene == 'VKORC1':
            warfarin_guidance = self.get_vkorc1_warfarin_guidance(sorted_variants)
        elif gene == 'IFNL3':
            hcv_guidance = self.get_ifnl3_hcv_response(sorted_variants)
        
        return NonStarAlleleFinding(
            gene=gene,
            variants=sorted_variants,
            phenotype_predictions={},
            drug_recommendations=drug_recommendations,
            warfarin_dosing_guidance=warfarin_guidance,
            hepatitis_response=hcv_guidance,
            clinical_significance=clinical_significance
        )
    
    def _determine_clinical_significance(self, gene: str, variants: List[VariantAnnotation]) -> str:
        """Determine overall clinical significance"""
        if not variants:
            return 'LOW'
        
        # Check if any HIGH impact variants
        high_impact = any(v.impact == 'HIGH' for v in variants)
        if high_impact:
            return 'HIGH'
        
        # Check MODERATE impact
        moderate_impact = any(v.impact == 'MODERATE' for v in variants)
        if moderate_impact:
            return 'MODERATE'
        
        return 'LOW'
    
    def _extract_drug_recommendations(self, 
                                     gene: str, 
                                     variants: List[VariantAnnotation]) -> Dict[str, List[Dict[str, Any]]]:
        """Extract drug-specific recommendations"""
        recommendations = {}
        
        if gene in self.DRUG_RECOMMENDATIONS:
            base_recs = self.DRUG_RECOMMENDATIONS[gene]
            
            for drug, rec in base_recs.items():
                recommendations[drug] = []
                
                if isinstance(rec, dict):
                    for condition, guidance in rec.items():
                        recommendations[drug].append({
                            'condition': condition,
                            'guidance': guidance,
                            'evidence_level': 'Moderate'
                        })
                elif isinstance(rec, str):
                    recommendations[drug].append({
                        'condition': 'General',
                        'guidance': rec,
                        'evidence_level': 'Moderate'
                    })
        
        return recommendations
    
    def format_variant_report(self, finding: NonStarAlleleFinding) -> str:
        """
        Format finding for clinical report
        
        Args:
            finding: NonStarAlleleFinding object
            
        Returns:
            Formatted text report
        """
        report_lines = []
        
        # Header
        report_lines.append("=" * 80)
        report_lines.append(f"NON-STAR ALLELE GENE ANALYSIS: {finding.gene}")
        report_lines.append("=" * 80)
        
        # Clinical significance
        report_lines.append(f"\nClinical Significance: {finding.clinical_significance}")
        
        # Detected variants
        report_lines.append("\nDETECTED VARIANTS:")
        report_lines.append("-" * 80)
        for i, var in enumerate(finding.variants, 1):
            report_lines.append(f"\n{i}. {var.chromosome}:{var.position} {var.ref_allele}→{var.alt_allele}")
            if var.rsid:
                report_lines.append(f"   rsID: {var.rsid}")
            report_lines.append(f"   Impact: {var.impact}")
            report_lines.append(f"   Consequence: {var.consequence}")
            if var.genotype:
                report_lines.append(f"   Genotype: {var.genotype}")
            if var.allele_frequency:
                report_lines.append(f"   Allele Frequency: {var.allele_frequency:.4f}")
        
        # Drug recommendations
        if finding.drug_recommendations:
            report_lines.append("\n\nDRUG RECOMMENDATIONS:")
            report_lines.append("-" * 80)
            for drug, recs in finding.drug_recommendations.items():
                report_lines.append(f"\n• {drug.upper()}:")
                for rec in recs:
                    report_lines.append(f"  - {rec['guidance']}")
                    report_lines.append(f"    Evidence Level: {rec['evidence_level']}")
        
        # Gene-specific guidance
        if finding.warfarin_dosing_guidance:
            report_lines.append(f"\n\n{finding.warfarin_dosing_guidance}")
        
        if finding.hepatitis_response:
            report_lines.append(f"\n\n{finding.hepatitis_response}")
        
        report_lines.append("\n" + "=" * 80)
        
        return "\n".join(report_lines)
    
    def to_dict(self, finding: NonStarAlleleFinding) -> Dict[str, Any]:
        """Convert finding to dictionary for JSON output"""
        return {
            'gene': finding.gene,
            'clinical_significance': finding.clinical_significance,
            'variants': [
                {
                    'chromosome': v.chromosome,
                    'position': v.position,
                    'ref_allele': v.ref_allele,
                    'alt_allele': v.alt_allele,
                    'rsid': v.rsid,
                    'impact': v.impact,
                    'consequence': v.consequence,
                    'genotype': v.genotype,
                    'allele_frequency': v.allele_frequency
                }
                for v in finding.variants
            ],
            'drug_recommendations': {
                drug: [
                    {
                        'condition': rec['condition'],
                        'guidance': rec['guidance'],
                        'evidence_level': rec['evidence_level']
                    }
                    for rec in recs
                ]
                for drug, recs in finding.drug_recommendations.items()
            },
            'warfarin_dosing_guidance': finding.warfarin_dosing_guidance,
            'hepatitis_c_response': finding.hepatitis_response
        }
