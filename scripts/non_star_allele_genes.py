#!/usr/bin/env python3
"""
Non-Star Allele Gene Handler
Handles genes that don't use star allele nomenclature (VKORC1, IFNL3, etc.)
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
    Handles analysis of genes without standard star allele nomenclature
    Primary targets: VKORC1, IFNL3, etc.
    """
    
    # Strictly genes that DO NOT use star alleles. 
    # Removed CYP4F2, ABCG2, NUDT15, SLCO1B1, G6PD as they do map to star variants in standard panels.
    NON_STAR_ALLELE_GENES = {
        'VKORC1', 'IFNL3', 'F2', 'F5', 'MTHFR', 'PAH'
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
        '-1639G>A': {
            'genotype_effects': {
                'G/G': {'description': '-1639G/G', 'warfarin_sensitivity': 'lower', 'dose_prediction': 'higher'},
                'G/A': {'description': '-1639G/A', 'warfarin_sensitivity': 'intermediate', 'dose_prediction': 'intermediate'},
                'A/A': {'description': '-1639A/A', 'warfarin_sensitivity': 'high', 'dose_prediction': 'lower'}
            },
            'gene_description': 'Vitamin K epoxide reductase complex subunit 1',
            'clinical_relevance': 'Major warfarin dosing determinant',
            'reference_url': 'https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/'
        }
    }
    
    # IFNL3 (IL28B) hepatitis C response
    IFNL3_VARIANTS = {
        '-1595G>A': {
            'genotype_effects': {
                'C/C': {'description': 'IL28B CC', 'response_rate': '70-80%', 'recommendation': 'Standard therapy'},
                'C/T': {'description': 'IL28B CT', 'response_rate': '40-50%', 'recommendation': 'Consider extended therapy'},
                'T/T': {'description': 'IL28B TT', 'response_rate': '25-30%', 'recommendation': 'Extended therapy or DAA'}
            },
            'gene_description': 'Interferon-lambda 3',
            'clinical_relevance': 'Hepatitis C spontaneous clearance and treatment response',
            'reference_url': 'https://www.ncbi.nlm.nih.gov/grc/human/genes/282617'
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
        }
    }
    
    def __init__(self):
        """Initialize the handler"""
        self.pharmgkb_annotations: Dict[str, List[PharmGKBVariantAnnotation]] = {}
        
    def is_non_star_allele_gene(self, gene: str) -> bool:
        """Check if gene uses non-star-allele nomenclature"""
        return gene in self.NON_STAR_ALLELE_GENES
    
    def get_vkorc1_warfarin_guidance(self, mapped_alleles: List[str]) -> str:
        """Generate warfarin dosing guidance for VKORC1 using mapped variants"""
        guidance_lines = ["WARFARIN DOSING GUIDANCE (VKORC1):\n"]
        
        for allele in mapped_alleles:
            clean_allele = allele.replace("_hom", "")
            if clean_allele in self.VKORC1_VARIANTS:
                info = self.VKORC1_VARIANTS[clean_allele]
                genotype_key = 'A/A' if allele.endswith('_hom') else 'G/A'
                
                if genotype_key in info['genotype_effects']:
                    effect = info['genotype_effects'][genotype_key]
                    guidance_lines.append(f"  Detected {clean_allele} ({effect['description']}):")
                    guidance_lines.append(f"    - Warfarin Sensitivity: {effect['warfarin_sensitivity'].upper()}")
                    guidance_lines.append(f"    - Dose Prediction: {effect['dose_prediction'].upper()} doses")
                    guidance_lines.append(f"    - Expected INR: Monitor closely (target 2-3)")
                    return "\n".join(guidance_lines)

        guidance_lines.append("  Wild-Type / Normal - use standard warfarin dosing")
        return "\n".join(guidance_lines)
    
    def get_ifnl3_hcv_response(self, mapped_alleles: List[str]) -> str:
        """Generate HCV treatment response prediction for IFNL3"""
        guidance_lines = ["HEPATITIS C TREATMENT RESPONSE (IFNL3):\n"]
        
        for allele in mapped_alleles:
            clean_allele = allele.replace("_hom", "")
            if clean_allele in self.IFNL3_VARIANTS:
                info = self.IFNL3_VARIANTS[clean_allele]
                genotype_key = 'T/T' if allele.endswith('_hom') else 'C/T'
                
                if genotype_key in info['genotype_effects']:
                    effect = info['genotype_effects'][genotype_key]
                    guidance_lines.append(f"  Variant {clean_allele} - Genotype {genotype_key}:")
                    guidance_lines.append(f"    - Phenotype: {effect['description']}")
                    guidance_lines.append(f"    - HCV Response Rate: {effect['response_rate']}")
                    guidance_lines.append(f"    - Recommendation: {effect['recommendation']}")
                    return "\n".join(guidance_lines)
        
        guidance_lines.append("  No major IFNL3/IL28B variants detected (Standard CC Response)")
        return "\n".join(guidance_lines)
    
    def analyze_gene_variants(self, 
                             gene: str, 
                             variants: List[VariantAnnotation],
                             pharmgkb_annotations: Dict[str, List[str]] = None,
                             mapped_alleles: List[str] = None) -> NonStarAlleleFinding:
        """
        Comprehensive analysis of non-star-allele gene
        """
        sorted_variants = sorted(
            variants,
            key=lambda v: self.IMPACT_PRIORITY.get(v.impact, 99)
        )
        
        clinical_significance = self._determine_clinical_significance(gene, sorted_variants)
        drug_recommendations = self._extract_drug_recommendations(gene, sorted_variants)
        
        warfarin_guidance = None
        hcv_guidance = None
        
        mapped_alleles = mapped_alleles or []
        
        if gene == 'VKORC1':
            warfarin_guidance = self.get_vkorc1_warfarin_guidance(mapped_alleles)
        elif gene == 'IFNL3':
            hcv_guidance = self.get_ifnl3_hcv_response(mapped_alleles)
            
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
        
        if any(v.impact == 'HIGH' for v in variants):
            return 'HIGH'
        if any(v.impact == 'MODERATE' for v in variants):
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
        """Format finding for clinical report"""
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append(f"NON-STAR ALLELE GENE ANALYSIS: {finding.gene}")
        report_lines.append("=" * 80)
        
        report_lines.append(f"\nClinical Significance: {finding.clinical_significance}")
        
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
        
        if finding.drug_recommendations:
            report_lines.append("\n\nDRUG RECOMMENDATIONS:")
            report_lines.append("-" * 80)
            for drug, recs in finding.drug_recommendations.items():
                report_lines.append(f"\n• {drug.upper()}:")
                for rec in recs:
                    report_lines.append(f"  - {rec['guidance']}")
                    report_lines.append(f"    Evidence Level: {rec['evidence_level']}")
        
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