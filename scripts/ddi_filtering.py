#!/usr/bin/env python3
"""
DDI Filtering Module
Implements intelligent drug-drug interaction filtering for clinical relevance
"""

from typing import List, Dict, Any, Set
from enum import Enum


class DDISeverity(Enum):
    """DDI severity levels"""
    SEVERE = "SEVERE"
    MAJOR = "MAJOR"
    MODERATE = "MODERATE"
    MINOR = "MINOR"
    UNKNOWN = "UNKNOWN"
    
    def get_priority(self) -> int:
        """Get numerical priority (higher = more severe)"""
        priorities = {
            "SEVERE": 5,
            "MAJOR": 4,
            "MODERATE": 3,
            "MINOR": 2,
            "UNKNOWN": 1
        }
        return priorities.get(self.value, 0)


class DDIFilter:
    """Filters and prioritizes drug-drug interactions"""
    
    # Top 200 most commonly prescribed drugs (by frequency in US)
    COMMON_DRUGS = {
        # Cardiovascular
        'atorvastatin', 'lisinopril', 'amlodipine', 'metoprolol', 'losartan', 'simvastatin',
        'rosuvastatin', 'clopidogrel', 'carvedilol', 'valsartan', 'warfarin', 'apixaban',
        'rivaroxaban', 'diltiazem', 'verapamil', 'pravastatin', 'enalapril', 'ramipril',
        'furosemide', 'hydrochlorothiazide', 'spironolactone',
        
        # Psychiatric
        'sertraline', 'escitalopram', 'citalopram', 'fluoxetine', 'paroxetine', 'bupropion',
        'venlafaxine', 'duloxetine', 'amitriptyline', 'trazodone', 'mirtazapine', 'quetiapine',
        'aripiprazole', 'risperidone', 'olanzapine', 'alprazolam', 'lorazepam', 'clonazepam',
        'diazepam', 'zolpidem', 'eszopiclone', 'atomoxetine', 'methylphenidate',
        
        # Pain/Inflammation
        'ibuprofen', 'naproxen', 'celecoxib', 'tramadol', 'hydrocodone', 'oxycodone',
        'gabapentin', 'pregabalin', 'meloxicam', 'diclofenac', 'acetaminophen',
        
        # Gastrointestinal
        'omeprazole', 'esomeprazole', 'lansoprazole', 'pantoprazole', 'rabeprazole',
        'famotidine', 'ranitidine', 'ondansetron',
        
        # Diabetes  
        'metformin', 'insulin', 'glipizide', 'glyburide', 'sitagliptin', 'pioglitazone',
        'empagliflozin', 'liraglutide', 'dulaglutide',
        
        # Respiratory
        'albuterol', 'fluticasone', 'budesonide', 'montelukast', 'ipratropium', 'tiotropium',
        
        # Antibiotics/Antivirals
        'amoxicillin', 'azithromycin', 'ciprofloxacin', 'levofloxacin', 'doxycycline',
        'cephalexin', 'metronidazole', 'valacyclovir', 'acyclovir',
        
        # Thyroid
        'levothyroxine', 'methimazole',
        
        # Immunosuppressants
        'prednisone', 'methylprednisolone', 'tacrolimus', 'cyclosporine', 'azathioprine',
        
        # Anticoagulation
        'aspirin', 'enoxaparin', 'heparin',
        
        # Oncology (common)
        'tamoxifen', 'anastrozole', 'letrozole', 'methotrexate', 'fluorouracil',
        'cyclophosphamide', 'cisplatin', 'irinotecan',
        
        # Urological
        'tamsulosin', 'finasteride', 'sildenafil', 'tadalafil',
        
        # Antiepileptics
        'lamotrigine', 'levetiracetam', 'topiramate', 'carbamazepine', 'phenytoin',
        'valproic acid', 'oxcarbazepine',
        
        # GI disorders
        'dicyclomine', 'hyoscyamine',
        
        # Allergy
        'cetirizine', 'fexofenadine', 'loratadine', 'diphenhydramine',
        
        # Muscle relaxants
        'cyclobenzaprine', 'baclofen', 'tizanidine',
        
        # Ophthalmology
        'latanoprost', 'timolol',
        
        # Dermatology
        'tretinoin', 'isotretinoin', 'adapalene',
        
        # Supplements (commonly considered)
        'vitamin d', 'vitamin b12', 'folic acid', 'iron', 'calcium'
    }
    
    # Gene-specific substrate drugs (drugs metabolized by specific genes)
    GENE_SUBSTRATES = {
        'CYP2D6': {
            'codeine', 'tramadol', 'hydrocodone', 'oxycodone', 'metoprolol', 'carvedilol',
            'propranolol', 'amitriptyline', 'nortriptyline', 'desipramine', 'imipramine',
            'paroxetine', 'fluoxetine', 'venlafaxine', 'risperidone', 'aripiprazole',
            'haloperidol', 'atomoxetine', 'dextromethorphan', 'tamoxifen', 'ondansetron',
            'metoclopramide', 'doxepin', 'clomipramine', 'trimipramine', 'duloxetine'
        },
        'CYP2C19': {
            'clopidogrel', 'citalopram', 'escitalopram', 'amitriptyline', 'clomipramine',
            'imipramine', 'omeprazole', 'esomeprazole', 'lansoprazole', 'pantoprazole',
            'voriconazole', 'diazepam', 'clobazam', 'phenytoin', 'warfarin', 'sertraline',
            'moclobemide', 'carisoprodol', 'propranolol', 'cyclophosphamide'
        },
        'CYP2C9': {
            'warfarin', 'phenytoin', 'losartan', 'celecoxib', 'flurbiprofen', 'ibuprofen',
            'piroxicam', 'meloxicam', 'glipizide', 'glyburide', 'tolbutamide', 'diclofenac',
            'naproxen', 'indomethacin', 'valsartan', 'irbesartan', 'rosiglitazone',
            'pioglitazone', 'fluoxetine', 'fluvastatin', 'phenobarbital', 'valproic acid'
        },
        'CYP3A4': {
            'simvastatin', 'atorvastatin', 'lovastatin', 'tacrolimus', 'cyclosporine',
            'midazolam', 'triazolam', 'alprazolam', 'fentanyl', 'oxycodone', 'buspirone',
            'quetiapine', 'aripiprazole', 'venlafaxine', 'trazodone', 'erythromycin',
            'clarithromycin', 'sildenafil', 'tadalafil', 'amlodipine', 'diltiazem',
            'verapamil', 'nifedipine', 'apixaban', 'rivaroxaban', 'ticagrelor'
        },
        'CYP3A5': {
            'tacrolimus', 'cyclosporine', 'sirolimus', 'midazolam', 'nifedipine', 'vincristine'
        },
        'SLCO1B1': {
            'simvastatin', 'atorvastatin', 'pravastatin', 'rosuvastatin', 'pitavastatin',
            'lovastatin', 'fluvastatin', 'repaglinide', 'glyburide', 'rifampin', 'methotrexate'
        },
        'ABCG2': {
            'rosuvastatin', 'sulfasalazine', 'methotrexate', 'topotecan', 'imatinib',
            'sunitinib', 'gefitinib', 'allopurinol'
        },
        'DPYD': {
            'fluorouracil', 'capecitabine', 'tegafur'
        },
        'TPMT': {
            'mercaptopurine', 'azathioprine', 'thioguanine'
        },
        'NUDT15': {
            'mercaptopurine', 'azathioprine', 'thioguanine'
        },
        'UGT1A1': {
            'irinotecan', 'atazanavir', 'belinostat', 'nilotinib', 'etoposide'
        },
        'G6PD': {
            'rasburicase', 'primaquine', 'nitrofurantoin', 'dapsone', 'sulfamethoxazole',
            'aspirin', 'chloroquine', 'quinine'
        }
    }
    
    @classmethod
    def is_gene_specific(cls, drug: str, gene: str) -> bool:
        """Check if drug is a substrate for the gene
        
        Args:
            drug: Drug name
            gene: Gene name
        
        Returns:
            True if drug is a substrate for the gene
        """
        drug_lower = drug.lower().strip()
        gene_upper = gene.upper()
        
        if gene_upper not in cls.GENE_SUBSTRATES:
            return False
        
        gene_drugs = cls.GENE_SUBSTRATES[gene_upper]
        
        # Exact match
        if drug_lower in gene_drugs:
            return True
        
        # Partial match (for drug variations)
        for known_drug in gene_drugs:
            if known_drug in drug_lower or drug_lower in known_drug:
                return True
        
        return False
    
    @classmethod
    def is_commonly_prescribed(cls, drug: str) -> bool:
        """Check if drug is commonly prescribed
        
        Args:
            drug: Drug name
        
        Returns:
            True if commonly prescribed
        """
        drug_lower = drug.lower().strip()
        
        if drug_lower in cls.COMMON_DRUGS:
            return True
        
        # Partial match
        for common_drug in cls.COMMON_DRUGS:
            if common_drug in drug_lower or drug_lower in common_drug:
                return True
        
        return False
    
    @classmethod
    def infer_severity(cls, description: str) -> DDISeverity:
        """Infer interaction severity from description
        
        Args:
            description: Interaction description
        
        Returns:
            DDISeverity enum
        """
        if not description:
            return DDISeverity.UNKNOWN
        
        desc_lower = description.lower()
        
        # SEVERE keywords
        severe_keywords = [
            'contraindicated', 'avoid', 'do not', 'should not', 'must not',
            'fatal', 'life-threatening', 'severe', 'serious', 'significant increased risk',
            'black box', 'boxed warning', 'death', 'coma', 'seizure'
        ]
        if any(kw in desc_lower for kw in severe_keywords):
            return DDISeverity.SEVERE
        
        # MAJOR keywords
        major_keywords = [
            'increased risk', 'significantly', 'greatly', 'substantially',
            'major', 'important', 'closely monitor', 'dangerous', 'toxic',
            'toxicity', 'overdose', 'excessive', 'bleeding risk', 'hemorrhage'
        ]
        if any(kw in desc_lower for kw in major_keywords):
            return DDISeverity.MAJOR
        
        # MODERATE keywords
        moderate_keywords = [
            'increase', 'decrease', 'reduce', 'may affect', 'monitor',
            'caution', 'consider', 'adjust', 'potential', 'possible',
            'altered', 'enhanced', 'diminished'
        ]
        if any(kw in desc_lower for kw in moderate_keywords):
            return DDISeverity.MODERATE
        
        # DEFAULT to MINOR
        return DDISeverity.MINOR
    
    @classmethod
    def filter_ddis(cls, interactions: List[Dict[str, Any]], gene: str = None,
                   severity_threshold: DDISeverity = DDISeverity.MODERATE,
                   common_drugs_only: bool = True,
                   gene_specific_only: bool = True,
                   max_interactions: int = 200) -> List[Dict[str, Any]]:
        """Filter DDIs for clinical relevance
        
        Args:
            interactions: List of interaction dictionaries
            gene: Gene name for gene-specific filtering (optional)
            severity_threshold: Minimum severity to include
            common_drugs_only: Only include commonly prescribed drugs
            gene_specific_only: Only include gene-specific substrates
            max_interactions: Maximum number of interactions to return
        
        Returns:
            Filtered list of interactions
        """
        filtered = []
        
        for interaction in interactions:
            # Get severity
            sev_str = interaction.get('severity', 'UNKNOWN')
            try:
                severity = DDISeverity[sev_str]
            except KeyError:
                # Infer from description
                desc = interaction.get('description', '')
                severity = cls.infer_severity(desc)
                interaction['severity'] = severity.value
            
            # Check severity threshold
            if severity.get_priority() < severity_threshold.get_priority():
                continue
            
            # Get drug names
            drug_1 = (interaction.get('drug_1') or interaction.get('query_drug') or '').strip().lower()
            drug_2 = (interaction.get('drug_2') or interaction.get('interacting_drug_name') or '').strip().lower()
            
            # Gene-specific filtering
            if gene and gene_specific_only:
                if not (cls.is_gene_specific(drug_1, gene) or cls.is_gene_specific(drug_2, gene)):
                    continue
            
            # Common drugs filtering
            if common_drugs_only:
                if not (cls.is_commonly_prescribed(drug_1) or cls.is_commonly_prescribed(drug_2)):
                    continue
            
            filtered.append(interaction)
        
        # Sort by severity (highest first)
        filtered.sort(key=lambda x: DDISeverity[x.get('severity', 'UNKNOWN')].get_priority(), reverse=True)
        
        # Limit number
        return filtered[:max_interactions]
    
    @classmethod
    def deduplicate_ddis(cls, interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove duplicate interactions
        
        Args:
            interactions: List of interaction dictionaries
        
        Returns:
            Deduplicated list
        """
        seen = set()
        unique = []
        
        for interaction in interactions:
            drug_1 = (interaction.get('drug_1') or interaction.get('query_drug') or '').strip().lower()
            drug_2 = (interaction.get('drug_2') or interaction.get('interacting_drug_name') or '').strip().lower()
            desc = (interaction.get('description') or '').strip().lower()
            
            # Create signature
            signature = f"{min(drug_1, drug_2)}|{max(drug_1, drug_2)}|{desc[:100]}"
            
            if signature not in seen:
                seen.add(signature)
                unique.append(interaction)
        
        return unique


def summarize_ddi_statistics(interactions: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate statistics about DDI list
    
    Args:
        interactions: List of interaction dictionaries
    
    Returns:
        Statistics dictionary
    """
    if not interactions:
        return {
            'total': 0,
            'by_severity': {},
            'unique_drugs': 0
        }
    
    severity_counts = {}
    drugs = set()
    
    for interaction in interactions:
        # Count severity
        sev = interaction.get('severity', 'UNKNOWN')
        severity_counts[sev] = severity_counts.get(sev, 0) + 1
        
        # Collect unique drugs
        drug_1 = interaction.get('drug_1') or interaction.get('query_drug') or ''
        drug_2 = interaction.get('drug_2') or interaction.get('interacting_drug_name') or ''
        if drug_1:
            drugs.add(drug_1.lower())
        if drug_2:
            drugs.add(drug_2.lower())
    
    return {
        'total': len(interactions),
        'by_severity': severity_counts,
        'unique_drugs': len(drugs),
        'drugs': sorted(list(drugs))
    }


if __name__ == '__main__':
    # Test DDI filtering
    test_interactions = [
        {
            'drug_1': 'Warfarin',
            'drug_2': 'Aspirin',
            'description': 'The risk of bleeding can be significantly increased when Aspirin is combined with Warfarin.',
            'severity': 'MAJOR'
        },
        {
            'drug_1': 'Simvastatin',
            'drug_2': 'Amlodipine',
            'description': 'Amlodipine may increase serum levels of Simvastatin.',
            'severity': 'MODERATE'
        },
        {
            'drug_1': 'Some-Obscure-Drug-XYZ',
            'drug_2': 'Another-Rare-Drug',
            'description': 'Minor interaction',
            'severity': 'MINOR'
        }
    ]
    
    print("DDI Filtering Test")
    print("=" * 80)
    print(f"Original interactions: {len(test_interactions)}")
    
    filtered = DDIFilter.filter_ddis(
        test_interactions,
        gene='CYP2C9',
        severity_threshold=DDISeverity.MODERATE,
        common_drugs_only=True,
        gene_specific_only=False
    )
    
    print(f"Filtered interactions: {len(filtered)}")
    print("\nFiltered Results:")
    for interaction in filtered:
        print(f"  {interaction['drug_1']} + {interaction['drug_2']}: {interaction['severity']}")
    
    stats = summarize_ddi_statistics(filtered)
    print(f"\nStatistics: {stats}")
