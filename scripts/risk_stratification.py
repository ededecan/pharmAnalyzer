#!/usr/bin/env python3
"""
Risk Stratification Module
Implements color-coded risk priority system for pharmacogenomic reports
Based on commercial standards (Allelica, Genomind, MyMeds)
"""

from typing import Dict, Any, Optional, List
from enum import Enum


class RiskLevel(Enum):
    """Risk levels with color codes"""
    RED = "RED"              # Life-threatening, contraindicated
    ORANGE = "ORANGE"        # Abnormal, 1A evidence
    YELLOW = "YELLOW"        # Abnormal, 1B-3 evidence
    GREEN = "GREEN"          # Normal/standard dosing
    GRAY = "GRAY"            # Unknown/indeterminate
    
    def get_color_hex(self) -> str:
        """Get hexadecimal color code"""
        colors = {
            "RED": "#DC3545",
            "ORANGE": "#FD7E14",
            "YELLOW": "#FFC107",
            "GREEN": "#28A745",
            "GRAY": "#6C757D"
        }
        return colors.get(self.value, "#6C757D")
    
    def get_priority_score(self) -> int:
        """Get numerical priority (higher = more urgent)"""
        scores = {
            "RED": 5,
            "ORANGE": 4,
            "YELLOW": 3,
            "GREEN": 1,
            "GRAY": 2
        }
        return scores.get(self.value, 0)
    
    def get_description(self) -> str:
        """Get risk level description"""
        descriptions = {
            "RED": "Life-threatening toxicity risk - DO NOT USE",
            "ORANGE": "Abnormal allele with Level 1A evidence - High priority action required",
            "YELLOW": "Abnormal allele with Level 1B-3 evidence - Consider alternatives",
            "GREEN": "Normal allele - Standard dosing appropriate",
            "GRAY": "Unknown/Indeterminate - Clinical judgment required"
        }
        return descriptions.get(self.value, "Unknown risk level")


class RiskStratifier:
    """Calculates risk levels for gene-drug combinations"""
    
    # Life-threatening combinations (always RED)
    CRITICAL_COMBINATIONS = {
        ('HLA-B', '*15:02', 'carbamazepine'): 'Stevens-Johnson syndrome risk',
        ('HLA-B', '*15:02', 'oxcarbazepine'): 'Stevens-Johnson syndrome risk',
        ('HLA-B', '*15:02', 'phenytoin'): 'Stevens-Johnson syndrome risk',
        ('HLA-B', '*15:02', 'lamotrigine'): 'Stevens-Johnson syndrome risk',
        ('HLA-A', '*31:01', 'carbamazepine'): 'Drug reaction with eosinophilia',
        ('HLA-B', '*57:01', 'abacavir'): 'Hypersensitivity reaction',
        ('HLA-B', '*58:01', 'allopurinol'): 'Severe cutaneous adverse reactions',
        ('G6PD', 'Deficient', 'rasburicase'): 'Hemolytic anemia risk',
        ('G6PD', 'Deficient', 'primaquine'): 'Hemolytic anemia risk',
        ('DPYD', '*2A/*2A', 'fluorouracil'): 'Severe toxicity risk',
        ('DPYD', '*2A/*2A', 'capecitabine'): 'Severe toxicity risk',
        ('DPYD', '*13/*13', 'fluorouracil'): 'Severe toxicity risk',
        ('TPMT', 'PM', 'mercaptopurine'): 'Severe myelosuppression risk',
        ('TPMT', 'PM', 'azathioprine'): 'Severe myelosuppression risk',
        ('TPMT', 'PM', 'thioguanine'): 'Severe myelosuppression risk',
    }
    
    # Poor metabolizer drugs requiring dosage adjustment (ORANGE)
    PM_HIGH_PRIORITY_DRUGS = {
        'CYP2D6': ['codeine', 'tramadol', 'tamoxifen', 'atomoxetine', 'aripiprazole', 
                   'risperidone', 'haloperidol', 'metoprolol', 'carvedilol'],
        'CYP2C19': ['clopidogrel', 'citalopram', 'escitalopram', 'amitriptyline', 
                    'clomipramine', 'voriconazole'],
        'CYP2C9': ['warfarin', 'phenytoin', 'celecoxib', 'flurbiprofen'],
        'CYP3A5': ['tacrolimus'],
        'NUDT15': ['mercaptopurine', 'azathioprine', 'thioguanine'],
        'SLCO1B1': ['simvastatin', 'atorvastatin', 'pravastatin', 'rosuvastatin'],
        'UGT1A1': ['irinotecan', 'atazanavir']
    }
    
    @classmethod
    def calculate_risk(cls, gene: str, diplotype: str, phenotype: str, 
                      drug: str, evidence_level: str = None) -> RiskLevel:
        """Calculate risk level for gene-drug combination
        
        Args:
            gene: Gene name (e.g., 'CYP2D6', 'HLA-B')
            diplotype: Diplotype string (e.g., '*1/*2')
            phenotype: Phenotype/metabolizer status (e.g., 'Poor Metabolizer', 'PM')
            drug: Drug name
            evidence_level: PharmGKB evidence level ('1A', '1B', '2A', '2B', '3', '4')
        
        Returns:
            RiskLevel enum
        """
        gene = gene.upper()
        drug_lower = drug.lower()
        
        # Check for critical combinations (RED)
        for (crit_gene, crit_dip, crit_drug), _ in cls.CRITICAL_COMBINATIONS.items():
            if (gene == crit_gene.upper() and 
                cls._matches_diplotype(diplotype, crit_dip) and
                drug_lower == crit_drug.lower()):
                return RiskLevel.RED
        
        # Check phenotype-based risk
        phenotype_upper = phenotype.upper()
        
        # Poor Metabolizers with 1A evidence drugs (ORANGE)
        if 'POOR' in phenotype_upper or phenotype_upper == 'PM':
            if gene in cls.PM_HIGH_PRIORITY_DRUGS:
                if any(drug_lower == d.lower() for d in cls.PM_HIGH_PRIORITY_DRUGS[gene]):
                    if evidence_level == '1A':
                        return RiskLevel.ORANGE
                    elif evidence_level in ['1B', '2A']:
                        return RiskLevel.YELLOW
        
        # Ultrarapid Metabolizers with prodrugs (ORANGE for 1A)
        if 'ULTRARAPID' in phenotype_upper or 'ULTRA' in phenotype_upper or phenotype_upper == 'UM':
            prodrugs = ['codeine', 'tramadol', 'clopidogrel', 'tamoxifen']
            if any(drug_lower == d.lower() for d in prodrugs):
                if evidence_level == '1A':
                    return RiskLevel.ORANGE
                elif evidence_level in ['1B', '2A']:
                    return RiskLevel.YELLOW
        
        # Intermediate Metabolizers (YELLOW for 1A/1B)
        if 'INTERMEDIATE' in phenotype_upper or phenotype_upper == 'IM':
            if evidence_level in ['1A', '1B']:
                return RiskLevel.YELLOW
        
        # Normal/Standard (GREEN)
        if 'NORMAL' in phenotype_upper or phenotype_upper in ['NM', 'EXTENSIVE']:
            return RiskLevel.GREEN
        
        # Decreased function (YELLOW)
        if 'DECREASED' in phenotype_upper:
            if evidence_level in ['1A', '1B', '2A']:
                return RiskLevel.YELLOW
        
        # Unknown/Indeterminate (GRAY)
        if 'UNKNOWN' in phenotype_upper or 'INDETERMINATE' in phenotype_upper:
            return RiskLevel.GRAY
        
        # Default to GRAY
        return RiskLevel.GRAY
    
    @staticmethod
    def _matches_diplotype(diplotype: str, pattern: str) -> bool:
        """Check if diplotype matches a pattern (e.g., contains specific allele)"""
        if not diplotype or not pattern:
            return False
        diplotype_parts = diplotype.replace('/', ' ').replace('|', ' ').split()
        return pattern in diplotype_parts or pattern in diplotype
    
    @classmethod
    def get_critical_warning(cls, gene: str, diplotype: str, phenotype: str, drug: str) -> Optional[str]:
        """Get critical warning message if applicable
        
        Returns:
            Warning message or None
        """
        gene_upper = gene.upper()
        drug_lower = drug.lower()
        
        for (crit_gene, crit_dip, crit_drug), reason in cls.CRITICAL_COMBINATIONS.items():
            if (gene_upper == crit_gene.upper() and 
                cls._matches_diplotype(diplotype, crit_dip) and
                drug_lower == crit_drug.lower()):
                return f"⚠️ DO NOT USE {drug.upper()}: {reason}"
        
        return None
    
    @classmethod
    def get_action_recommendation(cls, risk_level: RiskLevel, gene: str, 
                                  phenotype: str, drug: str) -> str:
        """Get actionable recommendation based on risk level
        
        Returns:
            Clinical recommendation string
        """
        phenotype_upper = phenotype.upper()
        
        if risk_level == RiskLevel.RED:
            return f"DO NOT USE {drug}. Select alternative medication immediately."
        
        elif risk_level == RiskLevel.ORANGE:
            if 'POOR' in phenotype_upper or phenotype_upper == 'PM':
                if drug.lower() in ['clopidogrel']:
                    return f"Use alternative antiplatelet agent (prasugrel, ticagrelor). {drug} will have reduced efficacy."
                elif drug.lower() in ['codeine', 'tramadol']:
                    return f"Use alternative analgesic (morphine, hydromorphone, oxycodone). {drug} will have reduced efficacy."
                else:
                    return f"Consider 50% dose reduction or alternative medication. Monitor closely for adverse effects."
            elif 'ULTRA' in phenotype_upper or phenotype_upper == 'UM':
                return f"Consider dose increase or alternative medication. Monitor for reduced efficacy."
            return "Consult pharmacogenomics specialist or clinical pharmacist for dosing guidance."
        
        elif risk_level == RiskLevel.YELLOW:
            if 'INTERMEDIATE' in phenotype_upper or phenotype_upper == 'IM':
                return f"Consider 25-50% dose adjustment. Monitor therapeutic response and side effects."
            return "Use with caution. Monitor therapeutic drug levels if available."
        
        elif risk_level == RiskLevel.GREEN:
            return "Standard dosing appropriate. Follow standard prescribing guidelines."
        
        else:  # GRAY
            return "Limited evidence available. Use clinical judgment and monitor closely."


def prioritize_interactions(interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sort interactions by risk priority (highest risk first)
    
    Args:
        interactions: List of interaction dictionaries with 'risk_level' key
    
    Returns:
        Sorted list of interactions
    """
    def get_priority_score(interaction):
        risk_level = interaction.get('risk_level')
        if isinstance(risk_level, str):
            try:
                risk_level = RiskLevel[risk_level]
            except KeyError:
                risk_level = RiskLevel.GRAY
        elif not isinstance(risk_level, RiskLevel):
            risk_level = RiskLevel.GRAY
        return risk_level.get_priority_score()
    
    return sorted(interactions, key=get_priority_score, reverse=True)


if __name__ == '__main__':
    # Test risk stratification
    test_cases = [
        ('HLA-B', '*15:02/*01:01', 'Positive', 'carbamazepine', '1A'),
        ('CYP2C19', '*2/*2', 'Poor Metabolizer', 'clopidogrel', '1A'),
        ('CYP2D6', '*1/*2', 'Intermediate Metabolizer', 'codeine', '1A'),
        ('CYP2D6', '*1/*1', 'Normal Metabolizer', 'tramadol', '1A'),
        ('CYP2C9', '*1/*3', 'Intermediate Metabolizer', 'warfarin', '1B'),
    ]
    
    print("Risk Stratification Test Cases")
    print("=" * 80)
    
    for gene, diplotype, phenotype, drug, evidence in test_cases:
        risk = RiskStratifier.calculate_risk(gene, diplotype, phenotype, drug, evidence)
        action = RiskStratifier.get_action_recommendation(risk, gene, phenotype, drug)
        warning = RiskStratifier.get_critical_warning(gene, diplotype, phenotype, drug)
        
        print(f"\nGene: {gene} | Diplotype: {diplotype} | Drug: {drug}")
        print(f"Risk Level: {risk.value} ({risk.get_color_hex()})")
        print(f"Action: {action}")
        if warning:
            print(f"WARNING: {warning}")
        print("-" * 80)
