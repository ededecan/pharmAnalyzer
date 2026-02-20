#!/usr/bin/env python3
"""
Specific Dosing Guidance Module
Provides detailed dose adjustment percentages and specific dosing recommendations
based on pharmacogenomic phenotypes
"""

from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass
from enum import Enum


class DosageAdjustmentType(Enum):
    """Type of dosage adjustment"""
    PERCENTAGE = "Percentage reduction/increase"
    ABSOLUTE = "Absolute dose adjustment"
    RANGE = "Dose range"
    TITRATION = "Titration protocol"


@dataclass
class DrugDosing:
    """Detailed drug dosing information"""
    drug: str
    gene: str
    phenotype: str
    standard_dose: str
    adjusted_dose: str
    adjustment_type: DosageAdjustmentType
    adjustment_percentage: Optional[int]  # e.g., -25 for 25% reduction
    rationale: str
    monitoring_frequency: str
    warning_signs: List[str]
    contraindications: Optional[str]
    references: List[str]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'drug': self.drug,
            'gene': self.gene,
            'phenotype': self.phenotype,
            'standard_dose': self.standard_dose,
            'adjusted_dose': self.adjusted_dose,
            'adjustment_type': self.adjustment_type.value,
            'adjustment_percentage': self.adjustment_percentage,
            'rationale': self.rationale,
            'monitoring_frequency': self.monitoring_frequency,
            'warning_signs': self.warning_signs,
            'contraindications': self.contraindications
        }


class DosingGuidanceDatabase:
    """Comprehensive dosing guidance database"""
    
    DOSING_GUIDELINES = {
        # CYP2C19 - Citalopram/Escitalopram
        ('CYP2C19', 'PM', 'citalopram'): DrugDosing(
            drug='citalopram',
            gene='CYP2C19',
            phenotype='Poor Metabolizer',
            standard_dose='20-40 mg daily',
            adjusted_dose='10 mg daily initial, max 20 mg daily',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-50,
            rationale='PM cannot adequately metabolize citalopram; QT prolongation risk at higher doses',
            monitoring_frequency='Every 1 week x4, then bi-weekly; QTc monitoring at baseline and week 2',
            warning_signs=['Dizziness', 'Syncope', 'Palpitations', 'Severe nausea', 'Suicidal ideation'],
            contraindications='Do not exceed 20 mg daily; avoid with other QT-prolonging drugs',
            references=['PMID: 26059431', 'FDA Label 2011']
        ),
        
        ('CYP2C19', 'IM', 'citalopram'): DrugDosing(
            drug='citalopram',
            gene='CYP2C19',
            phenotype='Intermediate Metabolizer',
            standard_dose='20-40 mg daily',
            adjusted_dose='10-20 mg daily (standard lower end)',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-25,
            rationale='IM have reduced but not absent metabolism; use lower starting dose',
            monitoring_frequency='Every 1-2 weeks x4, then monthly',
            warning_signs=['Excessive sedation', 'Tremor', 'Agitation', 'Suicidal ideation'],
            contraindications='Max 30 mg daily recommended',
            references=['PMID: 26059431']
        ),
        
        ('CYP2C19', 'PM', 'escitalopram'): DrugDosing(
            drug='escitalopram',
            gene='CYP2C19',
            phenotype='Poor Metabolizer',
            standard_dose='10-20 mg daily',
            adjusted_dose='5 mg daily initial, max 10 mg daily',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-50,
            rationale='Active enantiomer; PM at high toxicity risk',
            monitoring_frequency='Every 1 week x4, then bi-weekly; QTc baseline and week 2',
            warning_signs=['QTc prolongation', 'Dizziness', 'Tremor', 'Weakness'],
            contraindications='Do not exceed 10 mg daily; monitor QTc if concurrent QT-prolonging drugs',
            references=['PMID: 26059431', 'FDA Label 2011']
        ),
        
        # CYP2D6 - Codeine (avoid in PM; adjust in UM)
        ('CYP2D6', 'UM', 'codeine'): DrugDosing(
            drug='codeine',
            gene='CYP2D6',
            phenotype='Ultrarapid Metabolizer',
            standard_dose='15-30 mg q4-6h PRN pain',
            adjusted_dose='AVOID or max 15 mg per dose, max 60 mg daily',
            adjustment_type=DosageAdjustmentType.ABSOLUTE,
            adjustment_percentage=None,
            rationale='Rapid conversion to morphine; risk of respiratory depression and death',
            monitoring_frequency='Avoid use; if necessary, close monitoring q4-6h for first 24h',
            warning_signs=['Respiratory depression', 'Drowsiness', 'Confusion', 'Loss of consciousness'],
            contraindications='AVOID codeine in UM; use non-opioid analgesics (acetaminophen, NSAIDs, tramadol)',
            references=['PMID: 19692698', 'FDA Label 2012']
        ),
        
        ('CYP2D6', 'PM', 'codeine'): DrugDosing(
            drug='codeine',
            gene='CYP2D6',
            phenotype='Poor Metabolizer',
            standard_dose='15-30 mg q4-6h PRN pain',
            adjusted_dose='AVOID - ineffective due to no morphine formation',
            adjustment_type=DosageAdjustmentType.ABSOLUTE,
            adjustment_percentage=None,
            rationale='Cannot convert codeine to active morphine; inadequate analgesia',
            monitoring_frequency='Do not use',
            warning_signs=['Inadequate pain relief'],
            contraindications='AVOID codeine in PM; alternative: morphine, hydrocodone, tramadol, NSAIDs',
            references=['PMID: 19692698', 'FDA Label 2007']
        ),
        
        ('CYP2D6', 'PM', 'tramadol'): DrugDosing(
            drug='tramadol',
            gene='CYP2D6',
            phenotype='Poor Metabolizer',
            standard_dose='50-100 mg q4-6h, max 400 mg/day',
            adjusted_dose='50 mg q6h, max 200 mg/day; consider alternative',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-50,
            rationale='PM have reduced tramadol metabolism; decreased analgesia and increased side effects',
            monitoring_frequency='Every 1-2 weeks for first month; assess pain control',
            warning_signs=['Inadequate analgesia', 'GI upset', 'CNS depression', 'Seizures (rare)'],
            contraindications='Use caution; risk of serotonin syndrome with SSRIs',
            references=['PMID: 28460141']
        ),
        
        # TPMT - Thiopurines
        ('TPMT', 'PM', 'mercaptopurine'): DrugDosing(
            drug='mercaptopurine',
            gene='TPMT',
            phenotype='Poor Metabolizer',
            standard_dose='1.5-2.5 mg/kg/day',
            adjusted_dose='0.15-0.25 mg/kg/day (5-10 fold reduction)',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-90,
            rationale='PM accumulate toxic 6-TG metabolites; high myelosuppression risk',
            monitoring_frequency='CBC weekly x4, then every 2 weeks',
            warning_signs=['Neutropenia', 'Thrombocytopenia', 'Anemia', 'Infections', 'Bleeding'],
            contraindications='Strict dose reduction required; consider alternative therapy (biologics)',
            references=['PMID: 18087212', 'PMID: 24074871']
        ),
        
        ('TPMT', 'IM', 'mercaptopurine'): DrugDosing(
            drug='mercaptopurine',
            gene='TPMT',
            phenotype='Intermediate Metabolizer',
            standard_dose='1.5-2.5 mg/kg/day',
            adjusted_dose='0.75-1.25 mg/kg/day (30-50% reduction)',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-40,
            rationale='IM have intermediate TG accumulation; moderate myelosuppression risk',
            monitoring_frequency='CBC weekly x2, then every 2-4 weeks',
            warning_signs=['Cytopenias', 'Hepatotoxicity signs', 'Rash'],
            contraindications='Monitor closely for marrow suppression',
            references=['PMID: 18087212']
        ),
        
        # UGT1A1 - Irinotecan
        ('UGT1A1', 'PM', 'irinotecan'): DrugDosing(
            drug='irinotecan',
            gene='UGT1A1',
            phenotype='Poor Metabolizer (TA7/TA7)',
            standard_dose='125-250 mg/m² IV q3 weeks',
            adjusted_dose='50-100 mg/m² IV q3 weeks (50% reduction)',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-50,
            rationale='PM have severely elevated SN-38 concentrations; high severe toxicity risk',
            monitoring_frequency='CBC, LFTs weekly x4; assess for diarrhea/neutropenia',
            warning_signs=['Severe neutropenia (<1000)', 'Severe diarrhea', 'Vomiting', 'Nausea'],
            contraindications='Consider alternative chemotherapy; if used, close supportive care needed',
            references=['PMID: 17053204', 'FDA Label 2013']
        ),
        
        ('UGT1A1', 'IM', 'irinotecan'): DrugDosing(
            drug='irinotecan',
            gene='UGT1A1',
            phenotype='Intermediate Metabolizer (TA6/TA7)',
            standard_dose='125-250 mg/m² IV q3 weeks',
            adjusted_dose='100-200 mg/m² IV q3 weeks (25% reduction)',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-25,
            rationale='IM have moderately elevated SN-38; increased severe toxicity risk',
            monitoring_frequency='CBC weekly x4, then per protocol; diarrhea monitoring',
            warning_signs=['Neutropenia', 'Diarrhea', 'Nausea', 'Vomiting'],
            contraindications='Monitor closely; may need further dose reduction if severe toxicity',
            references=['PMID: 17053204']
        ),
        
        # CYP2C9/VKORC1 - Warfarin
        ('CYP2C9', '*2/*3', 'warfarin'): DrugDosing(
            drug='warfarin',
            gene='CYP2C9',
            phenotype='Poor Metabolizer (*2/*3)',
            standard_dose='5-7 mg daily (target INR 2-3)',
            adjusted_dose='2-3 mg daily initial; target INR 2-3',
            adjustment_type=DosageAdjustmentType.PERCENTAGE,
            adjustment_percentage=-50,
            rationale='PM have severely reduced warfarin metabolism; high bleeding risk',
            monitoring_frequency='INR baseline, day 3-5, day 7, week 2-4, then weekly until stable',
            warning_signs=['Easy bruising', 'Bleeding gums', 'Blood in stool/urine', 'Excessive bruising'],
            contraindications='Use pharmacogenomic dosing calculator; frequent INR monitoring essential',
            references=['PMID: 27488176', 'CPIC Guidelines']
        ),
        
        # SLCO1B1 - Simvastatin
        ('SLCO1B1', '*5/*5', 'simvastatin'): DrugDosing(
            drug='simvastatin',
            gene='SLCO1B1',
            phenotype='Poor transporter activity (*5/*5)',
            standard_dose='20-80 mg daily',
            adjusted_dose='10-20 mg daily maximum',
            adjustment_type=DosageAdjustmentType.ABSOLUTE,
            adjustment_percentage=-75,
            rationale='*5/*5 carriers have dramatically reduced statin clearance; myopathy risk',
            monitoring_frequency='CK baseline and if muscle symptoms; lipid panel q4-6 weeks',
            warning_signs=['Muscle pain', 'Weakness', 'Dark urine', 'Elevated CK'],
            contraindications='Avoid simvastatin; consider pravastatin or rosuvastatin instead',
            references=['PMID: 18809606', 'PMID: 31594736']
        ),
        
        # CYP3A5 - Tacrolimus
        ('CYP3A5', '*1/*1', 'tacrolimus'): DrugDosing(
            drug='tacrolimus',
            gene='CYP3A5',
            phenotype='Expressor (*1/*1)',
            standard_dose='0.1-0.2 mg/kg/day divided BID',
            adjusted_dose='May require 1.5-2x standard dose; dose individual by TDM',
            adjustment_type=DosageAdjustmentType.RANGE,
            adjustment_percentage=None,
            rationale='CYP3A5 expressors have rapid tacrolimus metabolism; standard doses subtherapeutic',
            monitoring_frequency='Tacrolimus levels (trough target 8-12 ng/mL post-transplant), every 3-5 days until stable',
            warning_signs=['Graft rejection', 'Subtherapeutic levels', 'Elevated creatinine (graft dysfunction)'],
            contraindications='Use therapeutic drug monitoring (TDM) to guide dosing',
            references=['PMID: 19141122', 'CPIC Guidelines']
        ),
        
        ('CYP3A5', '*3/*3', 'tacrolimus'): DrugDosing(
            drug='tacrolimus',
            gene='CYP3A5',
            phenotype='Non-expressor (*3/*3)',
            standard_dose='0.1-0.2 mg/kg/day divided BID',
            adjusted_dose='Standard dosing; may need dose reduction',
            adjustment_type=DosageAdjustmentType.RANGE,
            adjustment_percentage=None,
            rationale='CYP3A5 non-expressors have slower tacrolimus metabolism; standard doses often adequate',
            monitoring_frequency='Tacrolimus levels (trough target 8-12 ng/mL), every 3-5 days until stable',
            warning_signs=['Nephrotoxicity', 'Supratherapeutic levels', 'Elevated creatinine'],
            contraindications='Use TDM; watch for toxicity from accumulation',
            references=['PMID: 19141122']
        ),
    }
    
    @classmethod
    def get_dosing_guidance(cls, drug: str, gene: str, phenotype: str) -> Optional[DrugDosing]:
        """Get dosing guidance for drug-gene-phenotype combination"""
        key = (gene.upper(), phenotype, drug.lower())
        
        # Exact match
        if key in cls.DOSING_GUIDELINES:
            return cls.DOSING_GUIDELINES[key]
        
        # Partial match on phenotype
        for guid_key, dosing in cls.DOSING_GUIDELINES.items():
            if (guid_key[0].upper() == gene.upper() and 
                guid_key[2].lower() == drug.lower()):
                if phenotype.upper() in guid_key[1].upper():
                    return dosing
        
        return None
    
    @classmethod
    def get_drugs_for_gene(cls, gene: str) -> List[str]:
        """Get all drugs with dosing guidance for a gene"""
        drugs = set()
        for key in cls.DOSING_GUIDELINES.keys():
            if key[0].upper() == gene.upper():
                drugs.add(key[2])
        return sorted(list(drugs))
    
    @classmethod
    def format_dosing_guidance(cls, dosing: DrugDosing) -> str:
        """Format dosing guidance for clinical report"""
        lines = []
        
        lines.append("=" * 80)
        lines.append(f"DOSING GUIDANCE: {dosing.drug.upper()}")
        lines.append("=" * 80)
        
        lines.append(f"\nGene: {dosing.gene}")
        lines.append(f"Phenotype: {dosing.phenotype}")
        
        lines.append(f"\nStandard Dose: {dosing.standard_dose}")
        lines.append(f"Adjusted Dose: {dosing.adjusted_dose}")
        
        if dosing.adjustment_percentage:
            direction = "reduction" if dosing.adjustment_percentage < 0 else "increase"
            lines.append(f"Adjustment: {abs(dosing.adjustment_percentage)}% {direction}")
        
        lines.append(f"\nRationale: {dosing.rationale}")
        lines.append(f"\nMonitoring Frequency: {dosing.monitoring_frequency}")
        
        if dosing.warning_signs:
            lines.append("\nWarning Signs to Monitor:")
            for sign in dosing.warning_signs:
                lines.append(f"  • {sign}")
        
        if dosing.contraindications:
            lines.append(f"\nContraindications/Precautions:\n  {dosing.contraindications}")
        
        lines.append("\n" + "=" * 80)
        
        return "\n".join(lines)
    
    @classmethod
    def calculate_dose(cls, standard_dose_mg: float, 
                      adjustment_percentage: Optional[int]) -> Optional[float]:
        """Calculate adjusted dose from percentage"""
        if adjustment_percentage is None:
            return None
        
        multiplier = 1 + (adjustment_percentage / 100)
        return round(standard_dose_mg * multiplier, 1)


class PolypharmacyDosingConsiderations:
    """Handle drug-drug interactions affecting dosing"""
    
    CYP_INHIBITOR_SUBSTRATES = {
        'CYP3A4': {
            'inhibitors': ['ketoconazole', 'ritonavir', 'grapefruit juice'],
            'substrates': ['simvastatin', 'tacrolimus', 'cyclosporine'],
            'interaction_strength': 'STRONG'
        },
        'CYP2D6': {
            'inhibitors': ['fluoxetine', 'paroxetine', 'bupropion'],
            'substrates': ['codeine', 'tramadol', 'metoprolol'],
            'interaction_strength': 'STRONG'
        },
        'CYP2C9': {
            'inhibitors': ['fluconazole', 'NSAIDs', 'sulfamethoxazole'],
            'substrates': ['warfarin', 'phenytoin'],
            'interaction_strength': 'STRONG'
        },
        'CYP2C19': {
            'inhibitors': ['cimetidine', 'omeprazole'],
            'substrates': ['citalopram', 'escitalopram', 'warfarin'],
            'interaction_strength': 'MODERATE'
        }
    }
    
    @staticmethod
    def check_interactions(primary_drug: str, 
                          concomitant_drugs: List[str],
                          gene: str) -> List[Dict[str, str]]:
        """Check for CYP interactions affecting dosing"""
        interactions = []
        
        for other_drug in concomitant_drugs:
            if PolypharmacyDosingConsiderations.is_cyp_inhibitor(other_drug):
                interactions.append({
                    'drug': other_drug,
                    'type': 'CYP inhibition',
                    'recommendation': f'May increase {primary_drug} levels; consider dose reduction'
                })
            elif PolypharmacyDosingConsiderations.is_cyp_inducer(other_drug):
                interactions.append({
                    'drug': other_drug,
                    'type': 'CYP induction',
                    'recommendation': f'May decrease {primary_drug} levels; consider dose increase'
                })
        
        return interactions
    
    @staticmethod
    def is_cyp_inhibitor(drug: str) -> bool:
        """Check if drug is a CYP inhibitor"""
        inhibitors = set()
        for cyp_info in PolypharmacyDosingConsiderations.CYP_INHIBITOR_SUBSTRATES.values():
            inhibitors.update(cyp_info['inhibitors'])
        return drug.lower() in [i.lower() for i in inhibitors]
    
    @staticmethod
    def is_cyp_inducer(drug: str) -> bool:
        """Check if drug is a CYP inducer"""
        inducers = ['rifampin', 'carbamazepine', 'phenytoin', 'phenobarbital', 'st. john\'s wort']
        return drug.lower() in inducers
