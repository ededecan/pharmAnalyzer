#!/usr/bin/env python3
"""
DPWG (Dutch Pharmacogenetics Working Group) Recommendations Module
Provides DPWG guidelines that complement CPIC recommendations
Focus on European standards and additional gene-drug pairs
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from enum import Enum


class RecommendationLevel(Enum):
    """DPWG recommendation strength levels"""
    A = "A - Important"  # Important change in effect
    B = "B - Possible"    # Possible effect
    C = "C - Advise"      # Advice due to implications


class ActionType(Enum):
    """Types of actions recommended"""
    AVOID = "AVOID"
    NO_ACTION = "No action"
    REDUCE_DOSE = "Reduce dose"
    INCREASE_DOSE = "Increase dose"
    MONITOR = "Monitor"
    OTHER = "Other"


@dataclass
class DPWGRecommendation:
    """DPWG guideline recommendation"""
    gene: str
    drug: str
    phenotype: str
    action: ActionType
    recommendation_level: RecommendationLevel
    description: str
    guidance: str
    references: List[str]
    evidence_url: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'gene': self.gene,
            'drug': self.drug,
            'phenotype': self.phenotype,
            'action': self.action.value,
            'recommendation_level': self.recommendation_level.value,
            'description': self.description,
            'guidance': self.guidance,
            'evidence_url': self.evidence_url
        }


class DPWGGuidelines:
    """DPWG guideline database and lookup"""
    
    RECOMMENDATIONS = {
        # CYP2C19 - Citalopram (European focus)
        ('CYP2C19', 'PM', 'citalopram'): DPWGRecommendation(
            gene='CYP2C19',
            drug='citalopram',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='CYP2C19 poor metabolizers have reduced metabolism of citalopram',
            guidance='Reduce initial dose by 50%. Maximum recommended daily dose: 20 mg/day due to QT prolongation risk.',
            references=['PMID: 26059431', 'PMID: 28460141'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL1000/variant/?variantAlleleId=PA144951674'
        ),
        
        # CYP2C19 - Escitalopram
        ('CYP2C19', 'PM', 'escitalopram'): DPWGRecommendation(
            gene='CYP2C19',
            drug='escitalopram',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='CYP2C19 poor metabolizers have reduced metabolism of escitalopram',
            guidance='Reduce initial dose to 5 mg/day. Maximum dose: 10 mg/day.',
            references=['PMID: 26059431'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL1083089/variant/'
        ),
        
        # CYP2D6 - Codeine (no effect in PM)
        ('CYP2D6', 'PM', 'codeine'): DPWGRecommendation(
            gene='CYP2D6',
            drug='codeine',
            phenotype='Poor Metabolizer',
            action=ActionType.AVOID,
            recommendation_level=RecommendationLevel.A,
            description='CYP2D6 poor metabolizers cannot convert codeine to morphine',
            guidance='Avoid codeine use. Morphine is inadequately formed in PM. Choose alternative analgesia.',
            references=['PMID: 19692698', 'PMID: 19667657'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL47/variant/'
        ),
        
        # CYP2D6 - Tramadol
        ('CYP2D6', 'PM', 'tramadol'): DPWGRecommendation(
            gene='CYP2D6',
            drug='tramadol',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.B,
            description='CYP2D6 poor metabolizers have reduced metabolism of tramadol',
            guidance='Reduce dose or consider alternative analgesic. Monitor for decreased analgesic effect.',
            references=['PMID: 28460141'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL6438/variant/'
        ),
        
        # CYP2D6 - Venlafaxine
        ('CYP2D6', 'PM', 'venlafaxine'): DPWGRecommendation(
            gene='CYP2D6',
            drug='venlafaxine',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.B,
            description='Poor metabolizers have increased venlafaxine concentrations',
            guidance='Reduce dose by 25%. Monitor for adverse effects (nausea, dizziness, headache).',
            references=['PMID: 20699393'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL1131/variant/'
        ),
        
        # TPMT - Azathioprine
        ('TPMT', 'PM', 'azathioprine'): DPWGRecommendation(
            gene='TPMT',
            drug='azathioprine',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='TPMT PM accumulate toxic metabolites (6-TG)',
            guidance='Reduce dose by 5-10%. Increase monitoring for bone marrow suppression.',
            references=['PMID: 18087212', 'PMID: 24074871'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL2064/variant/'
        ),
        
        # TPMT - Thiopurine
        ('TPMT', 'PM', 'thiopurine'): DPWGRecommendation(
            gene='TPMT',
            drug='thiopurine',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='TPMT PM accumulate toxic 6-TG metabolites',
            guidance='Reduce dose 5-10 fold. Increase monitoring for bone marrow suppression.',
            references=['PMID: 18087212', 'PMID: 24074871'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL77982/variant/'
        ),
        
        # HLA-B*5701 - Abacavir
        ('HLA-B', '*5701', 'abacavir'): DPWGRecommendation(
            gene='HLA-B',
            drug='abacavir',
            phenotype='HLA-B*5701 positive',
            action=ActionType.AVOID,
            recommendation_level=RecommendationLevel.A,
            description='HLA-B*5701 carriers at high risk of abacavir hypersensitivity reaction',
            guidance='Avoid abacavir. Risk of potentially fatal hypersensitivity reaction.',
            references=['PMID: 18227476', 'FDA Label'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL1644/variant/'
        ),
        
        # UGT1A1 - Irinotecan
        ('UGT1A1', 'IM', 'irinotecan'): DPWGRecommendation(
            gene='UGT1A1',
            drug='irinotecan',
            phenotype='Intermediate Metabolizer',
            action=ActionType.MONITOR,
            recommendation_level=RecommendationLevel.A,
            description='UGT1A1 IM have elevated SN-38 (active metabolite) concentrations',
            guidance='Reduce dose by 25%. Monitor carefully for severe neutropenia and diarrhea.',
            references=['PMID: 17053204'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL2127/variant/'
        ),
        
        ('UGT1A1', 'PM', 'irinotecan'): DPWGRecommendation(
            gene='UGT1A1',
            drug='irinotecan',
            phenotype='Poor Metabolizer',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='UGT1A1 PM have severely elevated toxic metabolite concentrations',
            guidance='Reduce dose by 50%. Increase monitoring, consider alternative drug.',
            references=['PMID: 17053204', 'FDA Label'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL2127/variant/'
        ),
        
        # NAT2 - Isoniazid
        ('NAT2', 'PM', 'isoniazid'): DPWGRecommendation(
            gene='NAT2',
            drug='isoniazid',
            phenotype='Poor Metabolizer',
            action=ActionType.MONITOR,
            recommendation_level=RecommendationLevel.B,
            description='NAT2 slow acetylators have increased isoniazid concentrations',
            guidance='Monitor for peripheral neuropathy and hepatotoxicity. Consider lower dose.',
            references=['PMID: 17505640'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL1030/variant/'
        ),
        
        # SLCO1B1 - Simvastatin
        ('SLCO1B1', '*5/*5', 'simvastatin'): DPWGRecommendation(
            gene='SLCO1B1',
            drug='simvastatin',
            phenotype='Poor transporter activity (*5/*5)',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='SLCO1B1 *5/*5 carriers have reduced statin clearance',
            guidance='Reduce simvastatin dose to maximum 20 mg/day. Higher risk of myopathy.',
            references=['PMID: 18809606', 'PMID: 31594736'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL386633/variant/'
        ),
        
        # CYP3A5 - Calcineurin inhibitors (tacrolimus)
        ('CYP3A5', '*1/*1', 'tacrolimus'): DPWGRecommendation(
            gene='CYP3A5',
            drug='tacrolimus',
            phenotype='Expressor (*1/*1)',
            action=ActionType.INCREASE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='CYP3A5 expressors have increased tacrolimus metabolism',
            guidance='May require higher doses to maintain therapeutic levels. TDM recommended.',
            references=['PMID: 19141122'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL375089/variant/'
        ),
        
        # CYP2C9/VKORC1 - Warfarin (European population focus)
        ('CYP2C9', '*2/*2', 'warfarin'): DPWGRecommendation(
            gene='CYP2C9',
            drug='warfarin',
            phenotype='Intermediate Metabolizer (*2/*2)',
            action=ActionType.REDUCE_DOSE,
            recommendation_level=RecommendationLevel.A,
            description='CYP2C9 *2 carriers have reduced warfarin metabolism',
            guidance='Use lower maintenance dose. Start with 2.5-3 mg daily. Frequent INR monitoring.',
            references=['PMID: 27488176'],
            evidence_url='https://www.pharmgkb.org/chemical/CHEMBL278/variant/'
        ),
    }
    
    @classmethod
    def get_recommendations(cls, gene: str, phenotype: str, drug: str) -> Optional[DPWGRecommendation]:
        """
        Get DPWG recommendation for gene-phenotype-drug combination
        
        Args:
            gene: Gene name (e.g., 'CYP2D6')
            phenotype: Metabolizer phenotype or genotype
            drug: Drug name
            
        Returns:
            DPWGRecommendation or None
        """
        # Try exact match first
        key = (gene, phenotype, drug)
        if key in cls.RECOMMENDATIONS:
            return cls.RECOMMENDATIONS[key]
        
        # Try case-insensitive match
        for rec_key, rec_value in cls.RECOMMENDATIONS.items():
            if (rec_key[0].upper() == gene.upper() and 
                rec_key[2].lower() == drug.lower()):
                if rec_key[1] in phenotype or phenotype in rec_key[1]:
                    return rec_value
        
        return None
    
    @classmethod
    def get_drugs_for_gene(cls, gene: str) -> List[str]:
        """Get all drugs with DPWG recommendations for a gene"""
        drugs = set()
        for key in cls.RECOMMENDATIONS.keys():
            if key[0].upper() == gene.upper():
                drugs.add(key[2])
        return sorted(list(drugs))
    
    @classmethod
    def get_all_genes_with_recommendations(cls) -> List[str]:
        """Get all genes with DPWG recommendations"""
        genes = set()
        for key in cls.RECOMMENDATIONS.keys():
            genes.add(key[0])
        return sorted(list(genes))
    
    @classmethod
    def search_recommendations(cls, search_term: str) -> List[DPWGRecommendation]:
        """Search for recommendations by gene or drug name"""
        results = []
        search_term = search_term.lower()
        
        for key, rec in cls.RECOMMENDATIONS.items():
            if (search_term in key[0].lower() or 
                search_term in key[2].lower() or
                search_term in rec.description.lower()):
                results.append(rec)
        
        return results


class DPWGRecommendationComparer:
    """Compare and integrate DPWG vs CPIC recommendations"""
    
    @staticmethod
    def compare_recommendations(gene: str, 
                                drug: str, 
                                phenotype: str,
                                dpwg_rec: Optional[DPWGRecommendation],
                                cpic_rec: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Compare DPWG and CPIC recommendations
        
        Args:
            gene: Gene name
            drug: Drug name
            phenotype: Metabolizer phenotype
            dpwg_rec: DPWG recommendation
            cpic_rec: CPIC recommendation (dict format)
            
        Returns:
            Comparison results
        """
        comparison = {
            'gene': gene,
            'drug': drug,
            'phenotype': phenotype,
            'dpwg': dpwg_rec.to_dict() if dpwg_rec else None,
            'cpic': cpic_rec,
            'consensus': None,
            'differences': [],
            'notes': []
        }
        
        if dpwg_rec and cpic_rec:
            # Check for agreement
            dpwg_action = dpwg_rec.action.value
            cpic_action = cpic_rec.get('recommendation', '')
            
            if dpwg_action.lower() == cpic_action.lower():
                comparison['consensus'] = 'AGREE'
                comparison['notes'].append('DPWG and CPIC recommendations align')
            else:
                comparison['consensus'] = 'DIFFER'
                comparison['differences'].append(f"DPWG action: {dpwg_action}")
                comparison['differences'].append(f"CPIC action: {cpic_action}")
                comparison['notes'].append('Recommendations differ - clinical judgment needed')
        
        elif dpwg_rec:
            comparison['consensus'] = 'DPWG ONLY'
            comparison['notes'].append('DPWG guideline available; no CPIC guideline found')
        
        elif cpic_rec:
            comparison['consensus'] = 'CPIC ONLY'
            comparison['notes'].append('CPIC guideline available; no DPWG guideline found')
        
        else:
            comparison['consensus'] = 'NO GUIDELINES'
            comparison['notes'].append('Neither DPWG nor CPIC guidelines found')
        
        return comparison
