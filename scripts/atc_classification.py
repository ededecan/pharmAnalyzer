#!/usr/bin/env python3
"""
ATC Classification Module
Organizes drugs by therapeutic area using Anatomical Therapeutic Chemical (ATC) classification
"""

from typing import Dict, List, Set, Optional, Tuple
from enum import Enum


class TherapeuticArea(Enum):
    """Major therapeutic areas"""
    CARDIOLOGY = "Cardiology"
    PSYCHIATRY = "Psychiatry"
    ONCOLOGY = "Oncology"
    PAIN_MANAGEMENT = "Pain Management"
    NEUROLOGY = "Neurology"
    GASTROENTEROLOGY = "Gastroenterology"
    HEMATOLOGY = "Hematology"
    IMMUNOLOGY = "Immunology"
    INFECTIOUS_DISEASE = "Infectious Disease"
    ENDOCRINOLOGY = "Endocrinology"
    PULMONOLOGY = "Pulmonology"
    RHEUMATOLOGY = "Rheumatology"
    DERMATOLOGY = "Dermatology"
    UROLOGY = "Urology"
    OTHER = "Other"


class ATCClassifier:
    """Classifies drugs into therapeutic areas"""
    
    # Major drug categories by therapeutic area
    DRUG_CATEGORIES = {
        TherapeuticArea.CARDIOLOGY: {
            'anticoagulants': ['warfarin', 'apixaban', 'rivaroxaban', 'dabigatran', 'edoxaban', 'acenocoumarol', 'phenprocoumon'],
            'antiplatelets': ['clopidogrel', 'prasugrel', 'ticagrelor', 'aspirin', 'dipyridamole'],
            'statins': ['simvastatin', 'atorvastatin', 'pravastatin', 'rosuvastatin', 'fluvastatin', 'lovastatin', 'pitavastatin'],
            'beta_blockers': ['metoprolol', 'propranolol', 'carvedilol', 'bisoprolol', 'atenolol', 'nebivolol', 'labetalol'],
            'ace_inhibitors': ['enalapril', 'lisinopril', 'ramipril', 'perindopril', 'captopril', 'fosinopril', 'quinapril', 'benazepril'],
            'arbs': ['losartan', 'valsartan', 'irbesartan', 'candesartan', 'telmisartan', 'olmesartan'],
            'calcium_channel_blockers': ['amlodipine', 'diltiazem', 'verapamil', 'nifedipine', 'felodipine'],
            'diuretics': ['furosemide', 'hydrochlorothiazide', 'spironolactone', 'chlorthalidone', 'indapamide'],
            'antiarrhythmics': ['amiodarone', 'flecainide', 'propafenone', 'sotalol', 'dronedarone']
        },
        TherapeuticArea.PSYCHIATRY: {
            'ssris': ['citalopram', 'escitalopram', 'fluoxetine', 'paroxetine', 'sertraline', 'fluvoxamine'],
            'snris': ['venlafaxine', 'duloxetine', 'desvenlafaxine', 'levomilnacipran'],
            'tcas': ['amitriptyline', 'nortriptyline', 'imipramine', 'desipramine', 'clomipramine', 'doxepin', 'trimipramine'],
            'antipsychotics': ['risperidone', 'aripiprazole', 'quetiapine', 'olanzapine', 'clozapine', 'haloperidol', 'ziprasidone', 
                             'paliperidone', 'iloperidone', 'lurasidone', 'asenapine', 'brexpiprazole', 'cariprazine'],
            'mood_stabilizers': ['lithium', 'valproic acid', 'carbamazepine', 'lamotrigine', 'oxcarbazepine'],
            'anxiolytics': ['diazepam', 'lorazepam', 'alprazolam', 'clonazepam', 'buspirone', 'clobazam'],
            'stimulants': ['atomoxetine', 'methylphenidate', 'amphetamine', 'dexmethylphenidate', 'lisdexamfetamine'],
            'other_antidepressants': ['bupropion', 'mirtazapine', 'trazodone', 'vortioxetine', 'vilazodone']
        },
        TherapeuticArea.ONCOLOGY: {
            'chemotherapy': ['fluorouracil', 'capecitabine', 'irinotecan', 'oxaliplatin', 'cisplatin', 'carboplatin', 
                           'cyclophosphamide', 'methotrexate', 'doxorubicin', 'paclitaxel', 'docetaxel'],
            'targeted_therapy': ['imatinib', 'gefitinib', 'erlotinib', 'sorafenib', 'sunitinib', 'axitinib', 'pazopanib',
                               'vemurafenib', 'dabrafenib', 'trametinib', 'ibrutinib', 'venetoclax'],
            'immunotherapy': ['rituximab', 'trastuzumab', 'pembrolizumab', 'nivolumab', 'ipilimumab'],
            'hormone_therapy': ['tamoxifen', 'anastrozole', 'letrozole', 'exemestane', 'fulvestrant'],
            'antimetabolites': ['mercaptopurine', 'azathioprine', 'thioguanine', 'cladribine', 'fludarabine'],
            'topoisomerase_inhibitors': ['topotecan', 'etoposide']
        },
        TherapeuticArea.PAIN_MANAGEMENT: {
            'opioids': ['codeine', 'tramadol', 'hydrocodone', 'oxycodone', 'morphine', 'hydromorphone', 
                       'fentanyl', 'methadone', 'oxymorphone', 'tapentadol', 'buprenorphine'],
            'nsaids': ['ibuprofen', 'naproxen', 'celecoxib', 'diclofenac', 'indomethacin', 'ketorolac', 
                      'meloxicam', 'piroxicam', 'flurbiprofen', 'lornoxicam', 'tenoxicam'],
            'other_analgesics': ['acetaminophen', 'paracetamol', 'gabapentin', 'pregabalin']
        },
        TherapeuticArea.NEUROLOGY: {
            'antiepileptics': ['phenytoin', 'carbamazepine', 'valproic acid', 'lamotrigine', 'levetiracetam', 
                             'oxcarbazepine', 'topiramate', 'zonisamide', 'lacosamide', 'perampanel'],
            'parkinsons': ['levodopa', 'pramipexole', 'ropinirole', 'rasagiline', 'selegiline'],
            'migraine': ['sumatriptan', 'rizatriptan', 'eletriptan', 'topiramate', 'propranolol']
        },
        TherapeuticArea.GASTROENTEROLOGY: {
            'ppis': ['omeprazole', 'esomeprazole', 'lansoprazole', 'pantoprazole', 'rabeprazole', 'dexlansoprazole'],
            'h2_blockers': ['famotidine', 'ranitidine', 'cimetidine'],
            'antiemetics': ['ondansetron', 'granisetron', 'metoclopramide', 'prochlorperazine']
        },
        TherapeuticArea.HEMATOLOGY: {
            'anticoagulants': ['warfarin', 'heparin', 'enoxaparin', 'fondaparinux'],
            'antiplatelets': ['clopidogrel', 'ticagrelor', 'prasugrel'],
            'thrombolytics': ['alteplase', 'tenecteplase']
        },
        TherapeuticArea.IMMUNOLOGY: {
            'immunosuppressants': ['tacrolimus', 'cyclosporine', 'azathioprine', 'mycophenolate', 'sirolimus', 'everolimus'],
            'dmards': ['methotrexate', 'sulfasalazine', 'hydroxychloroquine', 'leflunomide'],
            'biologics': ['adalimumab', 'etanercept', 'infliximab', 'rituximab', 'tocilizumab']
        },
        TherapeuticArea.INFECTIOUS_DISEASE: {
            'antivirals': ['efavirenz', 'nevirapine', 'abacavir', 'tenofovir', 'lamivudine', 'emtricitabine', 'dolutegravir', 'raltegravir'],
            'antibiotics': ['ceftriaxone', 'clarithromycin', 'amoxicillin', 'azithromycin', 'ciprofloxacin'],
            'antifungals': ['voriconazole', 'fluconazole', 'itraconazole', 'posaconazole', 'isavuconazole'],
            'antimalarials': ['artemether', 'lumefantrine', 'primaquine', 'chloroquine']
        },
        TherapeuticArea.ENDOCRINOLOGY: {
            'diabetes': ['metformin', 'glyburide', 'glipizide', 'pioglitazone', 'rosiglitazone', 'nateglinide', 'repaglinide', 'sitagliptin'],
            'thyroid': ['levothyroxine', 'methimazole', 'propylthiouracil']
        },
        TherapeuticArea.PULMONOLOGY: {
            'bronchodilators': ['albuterol', 'salmeterol', 'formoterol', 'tiotropium'],
            'corticosteroids': ['budesonide', 'fluticasone', 'beclomethasone', 'prednisone']
        },
        TherapeuticArea.RHEUMATOLOGY: {
            'antigout': ['allopurinol', 'febuxostat', 'colchicine', 'probenecid'],
            'nsaids': ['ibuprofen', 'naproxen', 'celecoxib', 'indomethacin'],
            'dmards': ['methotrexate', 'sulfasalazine', 'hydroxychloroquine']
        },
        TherapeuticArea.DERMATOLOGY: {
            'retinoids': ['isotretinoin', 'tretinoin', 'adapalene'],
            'immunosuppressants': ['methotrexate', 'cyclosporine']
        }
    }
    
    # Build reverse lookup: drug -> (therapeutic_area, category)
    _DRUG_LOOKUP = {}
    
    @classmethod
    def _build_lookup(cls):
        """Build drug lookup dictionary"""
        if cls._DRUG_LOOKUP:
            return
        
        for area, categories in cls.DRUG_CATEGORIES.items():
            for category, drugs in categories.items():
                for drug in drugs:
                    cls._DRUG_LOOKUP[drug.lower()] = (area, category)
    
    @classmethod
    def classify_drug(cls, drug_name: str) -> Tuple[TherapeuticArea, str]:
        """Classify a drug into therapeutic area and category
        
        Args:
            drug_name: Drug name
        
        Returns:
            Tuple of (TherapeuticArea, category_name) or (OTHER, 'other')
        """
        cls._build_lookup()
        drug_lower = drug_name.lower().strip()
        
        # Direct lookup
        if drug_lower in cls._DRUG_LOOKUP:
            return cls._DRUG_LOOKUP[drug_lower]
        
        # Partial matching for drug variations
        for known_drug, (area, category) in cls._DRUG_LOOKUP.items():
            if known_drug in drug_lower or drug_lower in known_drug:
                return (area, category)
        
        return (TherapeuticArea.OTHER, 'other')
    
    @classmethod
    def get_therapeutic_area(cls, drug_name: str) -> TherapeuticArea:
        """Get therapeutic area for a drug
        
        Args:
            drug_name: Drug name
        
        Returns:
            TherapeuticArea enum
        """
        area, _ = cls.classify_drug(drug_name)
        return area
    
    @classmethod
    def organize_drugs_by_area(cls, drug_list: List[str]) -> Dict[TherapeuticArea, List[str]]:
        """Organize a list of drugs by therapeutic area
        
        Args:
            drug_list: List of drug names
        
        Returns:
            Dictionary mapping TherapeuticArea to list of drugs
        """
        organized = {area: [] for area in TherapeuticArea}
        
        for drug in drug_list:
            area = cls.get_therapeutic_area(drug)
            organized[area].append(drug)
        
        # Remove empty areas
        return {area: drugs for area, drugs in organized.items() if drugs}
    
    @classmethod
    def organize_drugs_detailed(cls, drug_list: List[str]) -> Dict[str, Dict[str, List[str]]]:
        """Organize drugs by therapeutic area and category
        
        Args:
            drug_list: List of drug names
        
        Returns:
            Nested dictionary: {area_name: {category: [drugs]}}
        """
        organized = {}
        
        for drug in drug_list:
            area, category = cls.classify_drug(drug)
            area_name = area.value
            
            if area_name not in organized:
                organized[area_name] = {}
            
            if category not in organized[area_name]:
                organized[area_name][category] = []
            
            organized[area_name][category].append(drug)
        
        return organized
    
    @classmethod
    def get_common_drugs_for_gene(cls, gene: str) -> Dict[TherapeuticArea, List[str]]:
        """Get commonly affected drugs for a specific gene, organized by therapeutic area
        
        Args:
            gene: Gene name (e.g., 'CYP2D6')
        
        Returns:
            Dictionary mapping TherapeuticArea to list of drugs
        """
        gene_upper = gene.upper()
        
        # Gene-specific drug lists (most commonly affected)
        gene_drug_map = {
            'CYP2D6': [
                'codeine', 'tramadol', 'hydrocodone', 'oxycodone',
                'metoprolol', 'carvedilol', 'propranolol',
                'amitriptyline', 'nortriptyline', 'desipramine', 'imipramine',
                'paroxetine', 'fluoxetine', 'venlafaxine',
                'risperidone', 'aripiprazole', 'haloperidol',
                'atomoxetine', 'tamoxifen'
            ],
            'CYP2C19': [
                'clopidogrel', 'citalopram', 'escitalopram',
                'amitriptyline', 'clomipramine', 'imipramine',
                'omeprazole', 'esomeprazole', 'lansoprazole', 'pantoprazole',
                'voriconazole', 'diazepam', 'clobazam'
            ],
            'CYP2C9': [
                'warfarin', 'phenytoin', 'losartan',
                'celecoxib', 'flurbiprofen', 'ibuprofen', 'piroxicam',
                'glipizide', 'glyburide', 'tolbutamide'
            ],
            'CYP3A4': [
                'simvastatin', 'atorvastatin', 'lovastatin',
                'tacrolimus', 'cyclosporine',
                'midazolam', 'triazolam', 'alprazolam',
                'fentanyl', 'oxycodone'
            ],
            'CYP3A5': [
                'tacrolimus', 'cyclosporine', 'sirolimus'
            ],
            'SLCO1B1': [
                'simvastatin', 'atorvastatin', 'pravastatin', 'rosuvastatin', 
                'pitavastatin'
            ],
            'DPYD': [
                'fluorouracil', 'capecitabine', 'tegafur'
            ],
            'TPMT': [
                'mercaptopurine', 'azathioprine', 'thioguanine'
            ],
            'NUDT15': [
                'mercaptopurine', 'azathioprine', 'thioguanine'
            ],
            'UGT1A1': [
                'irinotecan', 'atazanavir', 'belinostat'
            ],
            'G6PD': [
                'rasburicase', 'primaquine', 'nitrofurantoin', 'dapsone'
            ]
        }
        
        drugs = gene_drug_map.get(gene_upper, [])
        return cls.organize_drugs_by_area(drugs)


if __name__ == '__main__':
    # Test ATC classification
    test_drugs = [
        'warfarin', 'clopidogrel', 'citalopram', 'codeine', 
        'irinotecan', 'tacrolimus', 'omeprazole', 'simvastatin'
    ]
    
    print("ATC Classification Test")
    print("=" * 80)
    
    for drug in test_drugs:
        area, category = ATCClassifier.classify_drug(drug)
        print(f"{drug:20s} -> {area.value:20s} ({category})")
    
    print("\n" + "=" * 80)
    print("\nOrganized by Therapeutic Area:")
    organized = ATCClassifier.organize_drugs_by_area(test_drugs)
    for area, drugs in organized.items():
        print(f"\n{area.value}:")
        for drug in drugs:
            print(f"  - {drug}")
