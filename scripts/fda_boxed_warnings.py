#!/usr/bin/env python3
"""
FDA Pharmacogenomics Boxed Warnings Module
Integrates FDA PgX-specific black box warnings and critical safety alerts
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from enum import Enum


class WarningCategory(Enum):
    """Categories of FDA warnings"""
    SEVERE_TOXICITY = "Severe Toxicity"
    HYPSENSITIVITY = "Hypersensitivity Reaction"
    APLASTIC_ANEMIA = "Aplastic Anemia"
    EFFECTIVE_CONTRACEPTION = "Teratogenicity"
    CARDIOVASCULAR = "Cardiovascular Risk"
    PSYCHIATRIC = "Psychiatric Risk"
    HEPATOTOXICITY = "Hepatotoxicity"
    NEPHROTOXICITY = "Nephrotoxicity"
    OTHER = "Other"


@dataclass
class FDABoxedWarning:
    """FDA pharmacogenomics boxed warning"""
    drug: str
    gene: str
    phenotype: str
    warning_category: WarningCategory
    warning_text: str
    clinical_implications: str
    recommended_action: str
    affected_populations_percentage: Optional[str]
    fda_approval_date: Optional[str]
    label_section: str
    reference_url: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'drug': self.drug,
            'gene': self.gene,
            'phenotype': self.phenotype,
            'warning_category': self.warning_category.value,
            'warning_text': self.warning_text,
            'clinical_implications': self.clinical_implications,
            'recommended_action': self.recommended_action,
            'affected_populations_percentage': self.affected_populations_percentage,
            'fda_approval_date': self.fda_approval_date,
            'label_section': self.label_section,
            'reference_url': self.reference_url
        }


class FDABoxedWarningDatabase:
    """FDA pharmacogenomics boxed warnings database"""
    
    WARNINGS = {
        # Abacavir - HLA-B*5701 hypersensitivity
        ('abacavir', 'HLA-B', '*5701'): FDABoxedWarning(
            drug='abacavir',
            gene='HLA-B',
            phenotype='HLA-B*5701 positive',
            warning_category=WarningCategory.HYPSENSITIVITY,
            warning_text='WARNING: HYPERSENSITIVITY REACTION\nAbacavir can cause a serious, potentially fatal hypersensitivity reaction.',
            clinical_implications='Patients with HLA-B*5701 allele are at significantly higher risk (approximately 5-8% of patients) of experiencing a hypersensitivity reaction.',
            recommended_action='HLA-B*5701 testing is REQUIRED before abacavir use. AVOID abacavir in HLA-B*5701 positive patients.',
            affected_populations_percentage='~5-8% of general population carry HLA-B*5701; varies by ethnicity',
            fda_approval_date='2008',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Carbamazepine - HLA-A*31:01, HLA-B*15:02
        ('carbamazepine', 'HLA-A', '*31:01'): FDABoxedWarning(
            drug='carbamazepine',
            gene='HLA-A',
            phenotype='HLA-A*31:01 positive',
            warning_category=WarningCategory.SEVERE_TOXICITY,
            warning_text='WARNING: SEVERE DERMATOLOGIC REACTIONS\nCarbamazepine can cause severe dermatologic reactions including Stevens-Johnson Syndrome (SJS) and Toxic Epidermal Necrolysis (TEN).',
            clinical_implications='HLA-A*31:01 carriers have dramatically increased risk of SJS/TEN, particularly in populations of European, Japanese, Korean, and Chinese ancestry.',
            recommended_action='Genetic testing for HLA-A*31:01 and HLA-B*15:02 is RECOMMENDED before initiating carbamazepine. Consider alternative anticonvulsants.',
            affected_populations_percentage='2-5% in European populations; 2-8% in East Asian populations',
            fda_approval_date='2007',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        ('carbamazepine', 'HLA-B', '*15:02'): FDABoxedWarning(
            drug='carbamazepine',
            gene='HLA-B',
            phenotype='HLA-B*15:02 positive',
            warning_category=WarningCategory.SEVERE_TOXICITY,
            warning_text='WARNING: SEVERE DERMATOLOGIC REACTIONS\nCarbamazepine can cause severe dermatologic reactions including Stevens-Johnson Syndrome (SJS) and Toxic Epidermal Necrolysis (TEN).',
            clinical_implications='HLA-B*15:02 carriers have dramatically increased risk of SJS/TEN, particularly in Southeast Asian populations.',
            recommended_action='Genetic testing for HLA-B*15:02 is REQUIRED in patients with ancestry in Southeast Asia. AVOID carbamazepine in HLA-B*15:02 positive patients.',
            affected_populations_percentage='2-12% in Southeast Asian populations',
            fda_approval_date='2007',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Azathioprine - TPMT deficiency
        ('azathioprine', 'TPMT', 'PM'): FDABoxedWarning(
            drug='azathioprine',
            gene='TPMT',
            phenotype='Poor Metabolizer',
            warning_category=WarningCategory.APLASTIC_ANEMIA,
            warning_text='WARNING: MYELOSUPPRESSION\nAzathioprine and related thiopurines can cause severe bone marrow suppression, including aplastic anemia.',
            clinical_implications='TPMT poor metabolizers accumulate toxic 6-thioguanine metabolites, leading to early-onset severe myelosuppression at standard doses.',
            recommended_action='TPMT testing is REQUIRED before azathioprine initiation. Dose reduction of 5-10 fold needed in TPMT PM.',
            affected_populations_percentage='~0.3% have TPMT PM phenotype',
            fda_approval_date='2013',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Irinotecan - UGT1A1 deficiency
        ('irinotecan', 'UGT1A1', 'TA7/TA7'): FDABoxedWarning(
            drug='irinotecan',
            gene='UGT1A1',
            phenotype='Poor Metabolizer (TA7/TA7)',
            warning_category=WarningCategory.SEVERE_TOXICITY,
            warning_text='WARNING: SEVERE NEUTROPENIA AND DIARRHEA\nIrinotecan can cause severe neutropenia and diarrhea, particularly in patients with reduced UGT1A1 activity.',
            clinical_implications='UGT1A1 TA7/TA7 carriers (Gilbert Syndrome variant) have elevated SN-38 concentrations and increased risk of severe toxicity.',
            recommended_action='UGT1A1 testing is RECOMMENDED. Dose reduction of 25-50% is needed in UGT1A1 PM. Close monitoring required.',
            affected_populations_percentage='~10-14% have Gilbert Syndrome; 0.3-1% have severe UGT1A1 deficiency',
            fda_approval_date='2003',
            label_section='Dosage and Administration',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Warfarin (CYP2C9/VKORC1)
        ('warfarin', 'CYP2C9', '*2'): FDABoxedWarning(
            drug='warfarin',
            gene='CYP2C9',
            phenotype='Intermediate Metabolizer (*1/*2)',
            warning_category=WarningCategory.CARDIOVASCULAR,
            warning_text='WARNING: BLEEDING\nWarfarin can cause major or fatal bleeding complications.',
            clinical_implications='CYP2C9 *2 carriers have reduced warfarin metabolism, increasing bleeding risk at standard doses.',
            recommended_action='Pharmacogenomic-guided dosing is RECOMMENDED. Initial dose should be reduced. More frequent INR monitoring needed.',
            affected_populations_percentage='Variable by ethnicity (5-30% carry CYP2C9 loss-of-function alleles)',
            fda_approval_date='2007',
            label_section='Clinical Pharmacology',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        ('warfarin', 'VKORC1', '-1639G/A'): FDABoxedWarning(
            drug='warfarin',
            gene='VKORC1',
            phenotype='VKORC1 c.-1639G/A A allele carrier',
            warning_category=WarningCategory.CARDIOVASCULAR,
            warning_text='WARNING: BLEEDING\nWarfarin can cause major or fatal bleeding complications.',
            clinical_implications='VKORC1 -1639G/A heterozygotes and GG homozygotes require lower warfarin doses.',
            recommended_action='Pharmacogenomic-guided dosing is RECOMMENDED. Genotyping can reduce time to therapeutic INR and bleeding complications.',
            affected_populations_percentage='40-60% carry VKORC1 variant alleles',
            fda_approval_date='2007',
            label_section='Clinical Pharmacology',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Clopidogrel - CYP2C19 loss-of-function
        ('clopidogrel', 'CYP2C19', 'PM'): FDABoxedWarning(
            drug='clopidogrel',
            gene='CYP2C19',
            phenotype='Poor Metabolizer',
            warning_category=WarningCategory.CARDIOVASCULAR,
            warning_text='WARNING: REDUCED EFFECTIVENESS\nClopidogrel effectiveness may be reduced in patients with CYP2C19 loss-of-function variants.',
            clinical_implications='CYP2C19 PM cannot adequately convert clopidogrel to active metabolite, reducing antiplatelet effect and increasing thrombotic risk.',
            recommended_action='CYP2C19 testing is RECOMMENDED. Consider alternative antiplatelet agents (prasugrel, ticagrelor) in CYP2C19 PM patients.',
            affected_populations_percentage='2-5% have CYP2C19 PM phenotype; varies by ethnicity',
            fda_approval_date='2009',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Codeine - CYP2D6 deficiency/ultrarapid metabolism
        ('codeine', 'CYP2D6', 'PM'): FDABoxedWarning(
            drug='codeine',
            gene='CYP2D6',
            phenotype='Poor Metabolizer',
            warning_category=WarningCategory.OTHER,
            warning_text='WARNING: REDUCED ANALGESIC EFFECT\nCodeine may be ineffective in patients with CYP2D6 poor metabolizer phenotype.',
            clinical_implications='CYP2D6 PM cannot convert codeine to morphine, resulting in inadequate analgesia.',
            recommended_action='CYP2D6 testing is RECOMMENDED. AVOID codeine in CYP2D6 PM. Use alternative analgesics.',
            affected_populations_percentage='5-10% have CYP2D6 PM phenotype; varies by ethnicity',
            fda_approval_date='2007',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        ('codeine', 'CYP2D6', 'UM'): FDABoxedWarning(
            drug='codeine',
            gene='CYP2D6',
            phenotype='Ultrarapid Metabolizer',
            warning_category=WarningCategory.SEVERE_TOXICITY,
            warning_text='WARNING: RESPIRATORY DEPRESSION AND DEATH\nCodeine is contraindicated in patients with CYP2D6 ultrarapid metabolizer phenotype.',
            clinical_implications='CYP2D6 UM rapidly converts codeine to morphine, resulting in potentially fatal respiratory depression at standard doses.',
            recommended_action='CYP2D6 testing is RECOMMENDED. AVOID codeine in CYP2D6 UM. Use non-opioid analgesics or alternative opioids.',
            affected_populations_percentage='1-5% have CYP2D6 UM phenotype; higher in some ethnic groups',
            fda_approval_date='2012',
            label_section='Contraindications',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Simvastatin - SLCO1B1
        ('simvastatin', 'SLCO1B1', '*5/*5'): FDABoxedWarning(
            drug='simvastatin',
            gene='SLCO1B1',
            phenotype='Poor transporter activity (*5/*5)',
            warning_category=WarningCategory.HEPATOTOXICITY,
            warning_text='WARNING: MYOPATHY AND RHABDOMYOLYSIS\nSimvastatin can cause muscle pain, weakness, and rhabdomyolysis, particularly in SLCO1B1 poor transporters.',
            clinical_implications='SLCO1B1 *5/*5 carriers have significantly elevated simvastatin concentrations, increasing myopathy risk.',
            recommended_action='Consider SLCO1B1 testing. Maximum simvastatin dose is 20 mg/day in SLCO1B1 *5/*5. Monitor for muscle symptoms.',
            affected_populations_percentage='3-5% have SLCO1B1 *5/*5 genotype',
            fda_approval_date='2012',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
        
        # Flucloxacillin - HLA-B*5701
        ('flucloxacillin', 'HLA-B', '*5701'): FDABoxedWarning(
            drug='flucloxacillin',
            gene='HLA-B',
            phenotype='HLA-B*5701 positive',
            warning_category=WarningCategory.HEPATOTOXICITY,
            warning_text='WARNING: DRUG-INDUCED LIVER INJURY\nFlucloxacillin may cause severe hepatotoxicity in HLA-B*5701 carriers.',
            clinical_implications='HLA-B*5701 carriers have 100+ fold increase in flucloxacillin-induced liver injury risk.',
            recommended_action='HLA-B*5701 testing is RECOMMENDED, especially in at-risk populations. Consider alternative antibiotics.',
            affected_populations_percentage='~5-8% general population; higher in some ethnic groups',
            fda_approval_date='2009',
            label_section='Warnings and Precautions',
            reference_url='https://www.fda.gov/drugs/drug-approvals-and-databases/fda-approved-drugs'
        ),
    }
    
    @classmethod
    def get_warning(cls, drug: str, gene: str, phenotype: str) -> Optional[FDABoxedWarning]:
        """Get FDA warning for drug-gene-phenotype combination"""
        key = (drug.lower(), gene.upper(), phenotype)
        
        # Exact match
        if key in cls.WARNINGS:
            return cls.WARNINGS[key]
        
        # Try partial match
        for warn_key, warning in cls.WARNINGS.items():
            if (warn_key[0].lower() == drug.lower() and 
                warn_key[1].upper() == gene.upper()):
                if phenotype.upper() in warn_key[2].upper():
                    return warning
        
        return None
    
    @classmethod
    def get_warnings_for_drug(cls, drug: str) -> List[FDABoxedWarning]:
        """Get all warnings for a drug"""
        warnings = []
        for key, warning in cls.WARNINGS.items():
            if key[0].lower() == drug.lower():
                warnings.append(warning)
        return warnings
    
    @classmethod
    def get_warnings_for_gene(cls, gene: str) -> List[FDABoxedWarning]:
        """Get all warnings for a gene"""
        warnings = []
        for key, warning in cls.WARNINGS.items():
            if key[1].upper() == gene.upper():
                warnings.append(warning)
        return warnings
    
    @classmethod
    def search_warnings(cls, search_term: str) -> List[FDABoxedWarning]:
        """Search warnings by drug, gene, or category"""
        results = []
        search_term_lower = search_term.lower()
        
        for warning in cls.WARNINGS.values():
            if (search_term_lower in warning.drug.lower() or
                search_term_lower in warning.gene.lower() or
                search_term_lower in warning.warning_category.value.lower() or
                search_term_lower in warning.warning_text.lower()):
                
                if warning not in results:
                    results.append(warning)
        
        return results


class WarningRiskAssessment:
    """Assess patient-specific risk based on genotype and warnings"""
    
    @staticmethod
    def assess_drug_safety(drug: str, 
                          genotypes: Dict[str, str]) -> Dict[str, Any]:
        """
        Assess drug safety given patient genotypes
        
        Args:
            drug: Drug name
            genotypes: Dictionary of gene -> phenotype/genotype
            
        Returns:
            Risk assessment
        """
        assessment = {
            'drug': drug,
            'risk_level': 'GREEN',  # GREEN, YELLOW (caution), RED (avoid)
            'warnings': [],
            'recommendations': [],
            'genetic_risk_factors': []
        }
        
        for gene, phenotype in genotypes.items():
            warning = FDABoxedWarningDatabase.get_warning(drug, gene, phenotype)
            if warning:
                assessment['warnings'].append(warning.to_dict())
                assessment['genetic_risk_factors'].append({
                    'gene': gene,
                    'phenotype': phenotype,
                    'category': warning.warning_category.value
                })
                
                # Determine risk level
                if 'avoid' in warning.recommended_action.lower() or 'contraindicated' in warning.recommended_action.lower():
                    assessment['risk_level'] = 'RED'
                    assessment['recommendations'].append(f"AVOID {drug} in {phenotype} patients")
                elif 'caution' in warning.recommended_action.lower() or 'monitor' in warning.recommended_action.lower():
                    if assessment['risk_level'] != 'RED':
                        assessment['risk_level'] = 'YELLOW'
                    assessment['recommendations'].append(warning.recommended_action)
        
        return assessment
