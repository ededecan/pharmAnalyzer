#!/usr/bin/env python3
"""
CPIC Guidelines Integration Module
Provides Clinical Pharmacogenetics Implementation Consortium (CPIC) recommendations
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass


@dataclass
class CPICRecommendation:
    """CPIC guideline recommendation"""
    gene: str
    diplotype: str
    phenotype: str
    drug: str
    recommendation: str
    classification: str  # 'Strong', 'Moderate', 'Optional'
    alternative_drugs: List[str]
    dosing_adjustment: Optional[str]
    evidence_level: str  # 'A', 'B', 'C', 'D'
    guideline_url: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'gene': self.gene,
            'diplotype': self.diplotype,
            'phenotype': self.phenotype,
            'drug': self.drug,
            'recommendation': self.recommendation,
            'classification': self.classification,
            'alternative_drugs': self.alternative_drugs,
            'dosing_adjustment': self.dosing_adjustment,
            'evidence_level': self.evidence_level,
            'guideline_url': self.guideline_url
        }


class CPICGuidelines:
    """CPIC guideline database and lookup"""
    
    # Comprehensive CPIC guidelines (Level A evidence)
    GUIDELINES = {
        # CYP2C19 - Clopidogrel
        ('CYP2C19', 'UM', 'clopidogrel'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*17/*17 or *1/*17',
            phenotype='Ultrarapid Metabolizer',
            drug='clopidogrel',
            recommendation='Use standard dosing. Monitor for increased bleeding risk.',
            classification='Moderate',
            alternative_drugs=[],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/'
        ),
        ('CYP2C19', 'RM', 'clopidogrel'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*1/*17',
            phenotype='Rapid Metabolizer',
            drug='clopidogrel',
            recommendation='Use standard dosing.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/'
        ),
        ('CYP2C19', 'NM', 'clopidogrel'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*1/*1',
            phenotype='Normal Metabolizer',
            drug='clopidogrel',
            recommendation='Use standard dosing.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/'
        ),
        ('CYP2C19', 'IM', 'clopidogrel'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*1/*2, *1/*3, *2/*17',
            phenotype='Intermediate Metabolizer',
            drug='clopidogrel',
            recommendation='Alternative antiplatelet therapy (prasugrel, ticagrelor) is recommended if no contraindication.',
            classification='Strong',
            alternative_drugs=['prasugrel', 'ticagrelor'],
            dosing_adjustment='If clopidogrel used: 150 mg loading, 75 mg daily',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/'
        ),
        ('CYP2C19', 'PM', 'clopidogrel'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*2/*2, *2/*3, *3/*3',
            phenotype='Poor Metabolizer',
            drug='clopidogrel',
            recommendation='Alternative antiplatelet therapy (prasugrel, ticagrelor) is recommended.',
            classification='Strong',
            alternative_drugs=['prasugrel', 'ticagrelor'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/'
        ),
        
        # CYP2C19 - SSRIs (Citalopram/Escitalopram)
        ('CYP2C19', 'UM', 'citalopram'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*17/*17',
            phenotype='Ultrarapid Metabolizer',
            drug='citalopram',
            recommendation='Consider alternative SSRI not predominantly metabolized by CYP2C19. If citalopram used, monitor for decreased efficacy.',
            classification='Moderate',
            alternative_drugs=['sertraline', 'paroxetine', 'fluvoxamine'],
            dosing_adjustment='Consider 50% dose increase',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-selective-serotonin-reuptake-inhibitors-and-cyp2d6-and-cyp2c19/'
        ),
        ('CYP2C19', 'PM', 'citalopram'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*2/*2, *2/*3, *3/*3',
            phenotype='Poor Metabolizer',
            drug='citalopram',
            recommendation='Consider 50% dose reduction or alternative SSRI not predominantly metabolized by CYP2C19.',
            classification='Strong',
            alternative_drugs=['sertraline', 'paroxetine', 'fluvoxamine'],
            dosing_adjustment='Reduce dose by 50%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-selective-serotonin-reuptake-inhibitors-and-cyp2d6-and-cyp2c19/'
        ),
        ('CYP2C19', 'PM', 'escitalopram'): CPICRecommendation(
            gene='CYP2C19',
            diplotype='*2/*2, *2/*3, *3/*3',
            phenotype='Poor Metabolizer',
            drug='escitalopram',
            recommendation='Consider 50% dose reduction or alternative SSRI not predominantly metabolized by CYP2C19.',
            classification='Strong',
            alternative_drugs=['sertraline', 'paroxetine', 'fluvoxamine'],
            dosing_adjustment='Reduce dose by 50%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-selective-serotonin-reuptake-inhibitors-and-cyp2d6-and-cyp2c19/'
        ),
        
        # CYP2D6 - Codeine
        ('CYP2D6', 'UM', 'codeine'): CPICRecommendation(
            gene='CYP2D6',
            diplotype='Gene duplications',
            phenotype='Ultrarapid Metabolizer',
            drug='codeine',
            recommendation='Avoid codeine use due to potential for toxicity. Use alternative analgesic (morphine, hydromorphone, oxycodone).',
            classification='Strong',
            alternative_drugs=['morphine', 'hydromorphone', 'oxycodone', 'buprenorphine'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/'
        ),
        ('CYP2D6', 'PM', 'codeine'): CPICRecommendation(
            gene='CYP2D6',
            diplotype='*3/*4, *4/*4, *5/*5',
            phenotype='Poor Metabolizer',
            drug='codeine',
            recommendation='Avoid codeine use due to lack of efficacy. Use alternative analgesic (morphine, hydromorphone, oxycodone).',
            classification='Strong',
            alternative_drugs=['morphine', 'hydromorphone', 'oxycodone', 'buprenorphine'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/'
        ),
        ('CYP2D6', 'NM', 'codeine'): CPICRecommendation(
            gene='CYP2D6',
            diplotype='*1/*1',
            phenotype='Normal Metabolizer',
            drug='codeine',
            recommendation='Use standard dosing.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/'
        ),
        
        # CYP2D6 - Tramadol
        ('CYP2D6', 'UM', 'tramadol'): CPICRecommendation(
            gene='CYP2D6',
            diplotype='Gene duplications',
            phenotype='Ultrarapid Metabolizer',
            drug='tramadol',
            recommendation='Avoid tramadol use. Use alternative analgesic.',
            classification='Strong',
            alternative_drugs=['morphine', 'hydromorphone', 'oxycodone', 'buprenorphine'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/'
        ),
        ('CYP2D6', 'PM', 'tramadol'): CPICRecommendation(
            gene='CYP2D6',
            diplotype='*3/*4, *4/*4, *5/*5',
            phenotype='Poor Metabolizer',
            drug='tramadol',
            recommendation='Avoid tramadol use. Use alternative analgesic.',
            classification='Strong',
            alternative_drugs=['morphine', 'hydromorphone', 'oxycodone', 'buprenorphine'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/'
        ),
        
        # CYP2C9 - Warfarin
        ('CYP2C9', 'IM', 'warfarin'): CPICRecommendation(
            gene='CYP2C9',
            diplotype='*1/*2, *1/*3',
            phenotype='Intermediate Metabolizer',
            drug='warfarin',
            recommendation='Consider 25-50% dose reduction. Use pharmacogenetic-based dosing algorithms.',
            classification='Strong',
            alternative_drugs=['apixaban', 'rivaroxaban', 'dabigatran'],
            dosing_adjustment='Reduce initial dose by 25-50%, titrate to INR',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/'
        ),
        ('CYP2C9', 'PM', 'warfarin'): CPICRecommendation(
            gene='CYP2C9',
            diplotype='*2/*2, *2/*3, *3/*3',
            phenotype='Poor Metabolizer',
            drug='warfarin',
            recommendation='Consider 50-70% dose reduction. Use pharmacogenetic-based dosing algorithms. Alternative anticoagulant may be preferred.',
            classification='Strong',
            alternative_drugs=['apixaban', 'rivaroxaban', 'dabigatran'],
            dosing_adjustment='Reduce initial dose by 50-70%, titrate cautiously to INR',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/'
        ),
        
        # TPMT - Thiopurines
        ('TPMT', 'PM', 'mercaptopurine'): CPICRecommendation(
            gene='TPMT',
            diplotype='*2/*2, *2/*3A, *3A/*3A, *3A/*3C, *3C/*3C',
            phenotype='Poor Metabolizer',
            drug='mercaptopurine',
            recommendation='Reduce daily dose by 10-fold and reduce frequency to 3 times per week OR use alternative agent.',
            classification='Strong',
            alternative_drugs=['alternative antimetabolite therapy'],
            dosing_adjustment='Reduce dose by 90%, give 3x weekly',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/'
        ),
        ('TPMT', 'IM', 'mercaptopurine'): CPICRecommendation(
            gene='TPMT',
            diplotype='*1/*2, *1/*3A, *1/*3B, *1/*3C',
            phenotype='Intermediate Metabolizer',
            drug='mercaptopurine',
            recommendation='Reduce starting dose by 30-70% depending on indication. Monitor for hematologic toxicity.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment='Reduce dose by 30-70%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/'
        ),
        ('TPMT', 'PM', 'azathioprine'): CPICRecommendation(
            gene='TPMT',
            diplotype='*2/*2, *2/*3A, *3A/*3A',
            phenotype='Poor Metabolizer',
            drug='azathioprine',
            recommendation='Reduce daily dose by 10-fold and reduce frequency to 3 times per week OR use alternative agent.',
            classification='Strong',
            alternative_drugs=['alternative immunosuppressant therapy'],
            dosing_adjustment='Reduce dose by 90%, give 3x weekly',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/'
        ),
        ('TPMT', 'IM', 'azathioprine'): CPICRecommendation(
            gene='TPMT',
            diplotype='*1/*2, *1/*3A, *1/*3B, *1/*3C',
            phenotype='Intermediate Metabolizer',
            drug='azathioprine',
            recommendation='Reduce starting dose by 30-70%. Monitor for hematologic toxicity.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment='Reduce dose by 30-70%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/'
        ),
        
        # DPYD - Fluoropyrimidines
        ('DPYD', 'PM', 'fluorouracil'): CPICRecommendation(
            gene='DPYD',
            diplotype='*2A/*2A, *13/*13',
            phenotype='Poor Metabolizer',
            drug='fluorouracil',
            recommendation='Avoid fluorouracil-containing regimens. Use alternative drug or non-fluoropyrimidine regimen.',
            classification='Strong',
            alternative_drugs=['alternative chemotherapy regimen'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/'
        ),
        ('DPYD', 'IM', 'fluorouracil'): CPICRecommendation(
            gene='DPYD',
            diplotype='*1/*2A, *1/*13, c.1236G>A/c.1236G>A',
            phenotype='Intermediate Metabolizer',
            drug='fluorouracil',
            recommendation='Reduce starting dose by 50%, followed by dose titration based on toxicity.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment='Reduce starting dose by 50%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/'
        ),
        ('DPYD', 'PM', 'capecitabine'): CPICRecommendation(
            gene='DPYD',
            diplotype='*2A/*2A, *13/*13',
            phenotype='Poor Metabolizer',
            drug='capecitabine',
            recommendation='Avoid capecitabine. Use alternative drug.',
            classification='Strong',
            alternative_drugs=['alternative chemotherapy regimen'],
            dosing_adjustment=None,
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/'
        ),
        ('DPYD', 'IM', 'capecitabine'): CPICRecommendation(
            gene='DPYD',
            diplotype='*1/*2A, *1/*13',
            phenotype='Intermediate Metabolizer',
            drug='capecitabine',
            recommendation='Reduce starting dose by 50%, followed by dose titration based on toxicity.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment='Reduce starting dose by 50%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/'
        ),
        
        # UGT1A1 - Irinotecan
        ('UGT1A1', 'PM', 'irinotecan'): CPICRecommendation(
            gene='UGT1A1',
            diplotype='*28/*28, *80/*80',
            phenotype='Poor Metabolizer',
            drug='irinotecan',
            recommendation='Reduce starting dose by at least 30%. Consider alternative therapy if high-dose regimen planned.',
            classification='Strong',
            alternative_drugs=[],
            dosing_adjustment='Reduce starting dose by 30-50%',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-ugt1a1-and-atazanavir/'
        ),
        ('UGT1A1', 'IM', 'irinotecan'): CPICRecommendation(
            gene='UGT1A1',
            diplotype='*1/*28, *1/*80',
            phenotype='Intermediate Metabolizer',
            drug='irinotecan',
            recommendation='Consider reducing dose by 25% for medium/high-dose regimens. Monitor for toxicity.',
            classification='Moderate',
            alternative_drugs=[],
            dosing_adjustment='Consider 25% dose reduction',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-ugt1a1-and-atazanavir/'
        ),
        
        # SLCO1B1 - Simvastatin  
        ('SLCO1B1', 'Decreased Function', 'simvastatin'): CPICRecommendation(
            gene='SLCO1B1',
            diplotype='*5/*5, *5/*15, *15/*15',
            phenotype='Decreased Function',
            drug='simvastatin',
            recommendation='Use lower doses (≤20 mg daily) or alternative statin (pravastatin, rosuvastatin ≤20 mg daily, or fluvastatin).',
            classification='Strong',
            alternative_drugs=['pravastatin', 'rosuvastatin', 'fluvastatin', 'pitavastatin'],
            dosing_adjustment='Maximum 20 mg daily',
            evidence_level='A',
            guideline_url='https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/'
        ),
    }
    
    @classmethod
    def get_recommendation(cls, gene: str, phenotype: str, drug: str) -> Optional[CPICRecommendation]:
        """Get CPIC recommendation for gene-phenotype-drug combination
        
        Args:
            gene: Gene name (e.g., 'CYP2D6')
            phenotype: Metabolizer phenotype (e.g., 'Poor Metabolizer', 'PM')
            drug: Drug name
        
        Returns:
            CPICRecommendation object or None
        """
        # Normalize phenotype
        phenotype_map = {
            'ULTRARAPID METABOLIZER': 'UM',
            'ULTRA-RAPID METABOLIZER': 'UM',
            'ULTRARAPID': 'UM',
            'RAPID METABOLIZER': 'RM',
            'NORMAL METABOLIZER': 'NM',
            'INTERMEDIATE METABOLIZER': 'IM',
            'POOR METABOLIZER': 'PM',
            'DECREASED FUNCTION': 'Decreased Function',
            'DEFICIENT': 'PM'
        }
        
        phenotype_normalized = phenotype_map.get(phenotype.upper(), phenotype.upper())
        
        # Try exact match
        key = (gene.upper(), phenotype_normalized, drug.lower())
        if key in cls.GUIDELINES:
            return cls.GUIDELINES[key]
        
        # Try phenotype variants
        for pheno_key in ['UM', 'RM', 'NM', 'IM', 'PM', 'Decreased Function']:
            key = (gene.upper(), pheno_key, drug.lower())
            if key in cls.GUIDELINES:
                if phenotype_normalized in [pheno_key, pheno_key.replace(' ', '_')]:
                    return cls.GUIDELINES[key]
        
        return None
    
    @classmethod
    def get_all_drugs_for_gene(cls, gene: str) -> List[str]:
        """Get all drugs with CPIC guidelines for a gene
        
        Args:
            gene: Gene name
        
        Returns:
            List of drug names
        """
        drugs = set()
        gene_upper = gene.upper()
        
        for (g, _, drug) in cls.GUIDELINES.keys():
            if g == gene_upper:
                drugs.add(drug)
        
        return sorted(list(drugs))
    
    @classmethod
    def has_guideline(cls, gene: str, phenotype: str, drug: str) -> bool:
        """Check if CPIC guideline exists
        
        Args:
            gene: Gene name
            phenotype: Metabolizer phenotype
            drug: Drug name
        
        Returns:
            True if guideline exists
        """
        return cls.get_recommendation(gene, phenotype, drug) is not None


if __name__ == '__main__':
    # Test CPIC guidelines
    test_cases = [
        ('CYP2C19', 'Poor Metabolizer', 'clopidogrel'),
        ('CYP2D6', 'Poor Metabolizer', 'codeine'),
        ('TPMT', 'Poor Metabolizer', 'mercaptopurine'),
        ('DPYD', 'Intermediate Metabolizer', 'fluorouracil'),
    ]
    
    print("CPIC Guidelines Test")
    print("=" * 80)
    
    for gene, phenotype, drug in test_cases:
        rec = CPICGuidelines.get_recommendation(gene, phenotype, drug)
        if rec:
            print(f"\n{gene} {phenotype} + {drug.upper()}")
            print(f"Classification: {rec.classification}")
            print(f"Recommendation: {rec.recommendation}")
            if rec.alternative_drugs:
                print(f"Alternatives: {', '.join(rec.alternative_drugs)}")
            if rec.dosing_adjustment:
                print(f"Dosing: {rec.dosing_adjustment}")
            print("-" * 80)
        else:
            print(f"\nNo guideline found for {gene} {phenotype} + {drug}")
            print("-" * 80)
