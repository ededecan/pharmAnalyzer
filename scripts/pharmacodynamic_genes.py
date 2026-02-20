#!/usr/bin/env python3
"""
Pharmacodynamic Genes Module
Handles genes involved in drug efficacy, safety, and disease response
Not classic metabolizer genes, but important for personalized medicine
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from enum import Enum


class GeneCategory(Enum):
    """Categories of pharmacodynamic genes"""
    IMMUNE_RESPONSE = "Immune Response"
    HEMOSTASIS = "Hemostasis/Thrombosis"
    FOLATE_METABOLISM = "Folate Metabolism"
    ELECTROLYTE = "Electrolyte Balance"
    PSYCHIATRIC = "Psychiatric Response"
    METABOLIC = "Metabolic Risk"
    TUMOR_SUPPRESSION = "Tumor Suppression"
    HYPERSENSITIVITY = "Drug Hypersensitivity"


@dataclass
class PDGeneInfo:
    """Pharmacodynamic gene information"""
    gene: str
    category: GeneCategory
    full_name: str
    biological_function: str
    key_variants: List[str]
    clinical_relevance: str
    drug_list: List[str]
    phenotype_implications: Dict[str, str]
    

@dataclass
class PDGeneRecommendation:
    """Recommendation for pharmacodynamic gene variant"""
    gene: str
    variant_or_phenotype: str
    drug: str
    recommendation: str
    clinical_rationale: str
    monitoring_parameters: List[str]
    alternatives: List[str]


class PharmacodynamicGenes:
    """Database of important pharmacodynamic genes"""
    
    GENE_INFO = {
        # HLA genes - drug hypersensitivity
        'HLA-A': PDGeneInfo(
            gene='HLA-A',
            category=GeneCategory.HYPERSENSITIVITY,
            full_name='Human Leukocyte Antigen A',
            biological_function='Major histocompatibility complex antigen; immune recognition',
            key_variants=['*31:01', '*24:02', '*02:01'],
            clinical_relevance='Strong genetic predictor of drug hypersensitivity reactions (carbamazepine, allopurinol)',
            drug_list=['carbamazepine', 'allopurinol', 'nevirapine'],
            phenotype_implications={
                '*31:01': 'High risk carbamazepine-induced Stevens-Johnson Syndrome (SJS)/TEN',
                '*24:02': 'Increased risk of nevirapine hypersensitivity',
                'Negative': 'Reduced hypersensitivity risk'
            }
        ),
        
        'HLA-B': PDGeneInfo(
            gene='HLA-B',
            category=GeneCategory.HYPERSENSITIVITY,
            full_name='Human Leukocyte Antigen B',
            biological_function='Major histocompatibility complex antigen; immune recognition',
            key_variants=['*5701', '*15:02', '*57:01'],
            clinical_relevance='Critical predictor of abacavir hypersensitivity and carbamazepine SJS/TEN',
            drug_list=['abacavir', 'carbamazepine', 'allopurinol', 'flucloxacillin'],
            phenotype_implications={
                '*5701': 'CONTRAINDICATION: 5-8% risk of abacavir hypersensitivity reaction',
                '*15:02': 'High risk of carbamazepine SJS/TEN in Southeast Asian populations',
                '*57:01': 'Abacavir precaution required',
                'Negative': 'Reduced abacavir hypersensitivity risk'
            }
        ),
        
        'HLA-DRB1': PDGeneInfo(
            gene='HLA-DRB1',
            category=GeneCategory.HYPERSENSITIVITY,
            full_name='Human Leukocyte Antigen DRB1',
            biological_function='Major histocompatibility complex antigen; immune recognition',
            key_variants=['*01:01', '*03:01', '*07:01'],
            clinical_relevance='Associated with drug hypersensitivity reactions and multiple sclerosis drug response',
            drug_list=['NSAIDs', 'sulfasalazine', 'DMARDs'],
            phenotype_implications={
                'Positive': 'Monitor for hypersensitivity reactions to immunomodulators',
                'Negative': 'Lower hypersensitivity risk'
            }
        ),
        
        # Hemostasis genes - bleeding/clotting risk
        'F2': PDGeneInfo(
            gene='F2',
            category=GeneCategory.HEMOSTASIS,
            full_name='Coagulation Factor II (thrombin)',
            biological_function='Key serine protease in coagulation cascade; converts fibrinogen to fibrin',
            key_variants=['G20210A', '-1691G>A'],
            clinical_relevance='Associated with increased thrombosis risk; modifies warfarin response',
            drug_list=['warfarin', 'acetylsalicylic acid', 'anticoagulants'],
            phenotype_implications={
                '20210A/A': '2-3 fold increased thrombosis risk; lower warfarin requirement',
                '20210G/A': '1.5-2 fold increased thrombosis risk',
                '20210G/G': 'Baseline thrombosis risk'
            }
        ),
        
        'F5': PDGeneInfo(
            gene='F5',
            category=GeneCategory.HEMOSTASIS,
            full_name='Coagulation Factor V (Leiden)',
            biological_function='Cofactor in coagulation cascade; target of activated protein C',
            key_variants=['Leiden (1691G>A)', 'HR2', 'H1299R'],
            clinical_relevance='Leiden mutation causes 5-10 fold thrombosis risk; modifies warfarin requirement',
            drug_list=['warfarin', 'oral contraceptives', 'hormone replacement therapy'],
            phenotype_implications={
                'Leiden/Leiden': '5-10 fold thrombosis risk even without other factors',
                'Leiden/WT': '2-3 fold increased thrombosis risk',
                'WT/WT': 'Baseline thrombosis risk; caution with OCP/HRT'
            }
        ),
        
        # Folate metabolism
        'MTHFR': PDGeneInfo(
            gene='MTHFR',
            category=GeneCategory.FOLATE_METABOLISM,
            full_name='Methylenetetrahydrofolate Reductase',
            biological_function='Catalyzes conversion of 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate',
            key_variants=['C677T', 'A1298C'],
            clinical_relevance='C677T variant reduces enzyme activity; affects methotrexate response and side effects',
            drug_list=['methotrexate', 'fluorouracil', 'flunitrazepam', 'warfarin'],
            phenotype_implications={
                'TT': '35% enzyme activity; potential methotrexate toxicity risk; may need folic acid supplementation',
                'CT': '65% enzyme activity',
                'CC': 'Normal enzyme activity'
            }
        ),
        
        # Psychiatric response genes
        'BDNF': PDGeneInfo(
            gene='BDNF',
            category=GeneCategory.PSYCHIATRIC,
            full_name='Brain-Derived Neurotrophic Factor',
            biological_function='Neurotrophin supporting nervous system development and plasticity',
            key_variants=['Val66Met'],
            clinical_relevance='Val66Met affects antidepressant and cognitive behavioral therapy response',
            drug_list=['sertraline', 'paroxetine', 'fluoxetine', 'venlafaxine'],
            phenotype_implications={
                'Met/Met': 'Potentially slower antidepressant response; may benefit from augmentation',
                'Val/Met': 'Intermediate response',
                'Val/Val': 'Better response to SSRIs and psychotherapy'
            }
        ),
        
        'COMT': PDGeneInfo(
            gene='COMT',
            category=GeneCategory.PSYCHIATRIC,
            full_name='Catechol-O-Methyltransferase',
            biological_function='Degrades dopamine, norepinephrine, and epinephrine',
            key_variants=['Val158Met'],
            clinical_relevance='Met/Met carriers have reduced COMT activity; may affect ADHD drug response',
            drug_list=['methylphenidate', 'amphetamine', 'atomoxetine'],
            phenotype_implications={
                'Met/Met': 'Reduced dopamine catabolism; may need lower stimulant doses',
                'Val/Met': 'Intermediate catabolism',
                'Val/Val': 'Higher catabolism; may need higher stimulant doses'
            }
        ),
        
        # Metabolic risk genes
        'MC4R': PDGeneInfo(
            gene='MC4R',
            category=GeneCategory.METABOLIC,
            full_name='Melanocortin-4 Receptor',
            biological_function='G-protein coupled receptor regulating energy homeostasis',
            key_variants=['rs52820871', 'rs17782313'],
            clinical_relevance='Associated with antipsychotic-induced weight gain and metabolic syndrome',
            drug_list=['olanzapine', 'clozapine', 'risperidone', 'paliperidone'],
            phenotype_implications={
                'Risk allele': 'Higher risk of antipsychotic-induced obesity',
                'WT/WT': 'Lower metabolic risk with antipsychotics'
            }
        ),
        
        'HTR2C': PDGeneInfo(
            gene='HTR2C',
            category=GeneCategory.PSYCHIATRIC,
            full_name='5-Hydroxytryptamine Receptor 2C',
            biological_function='G-protein coupled serotonin receptor; involved in appetite regulation',
            key_variants=['Cys23Ser'],
            clinical_relevance='Associated with weight gain and metabolic effects of SSRIs and antipsychotics',
            drug_list=['fluoxetine', 'sertraline', 'antipsychotics'],
            phenotype_implications={
                'Ser/Ser': 'Risk allele for SSRI-induced weight gain',
                'Cys/Ser': 'Intermediate risk',
                'Cys/Cys': 'Lower weight gain risk'
            }
        ),
        
        # TP53 - tumor suppression (oncology)
        'TP53': PDGeneInfo(
            gene='TP53',
            category=GeneCategory.TUMOR_SUPPRESSION,
            full_name='Tumor Protein P53 (Li-Fraumeni)',
            biological_function='Transcription factor controlling cell cycle and apoptosis; tumor suppressor',
            key_variants=['R175H', 'R248Q', 'R273H'],
            clinical_relevance='TP53 germline mutations (Li-Fraumeni) require modified chemotherapy and surveillance',
            drug_list=['doxorubicin', 'topoisomerase inhibitors', 'radiation therapy'],
            phenotype_implications={
                'Mutation carrier': 'CRITICAL: Increased risk of treatment-related secondary malignancies; avoid certain agents',
                'WT': 'Standard oncology dosing'
            }
        ),
    }
    
    @classmethod
    def get_gene_info(cls, gene: str) -> Optional[PDGeneInfo]:
        """Get information about pharmacodynamic gene"""
        return cls.GENE_INFO.get(gene.upper())
    
    @classmethod
    def get_genes_in_category(cls, category: GeneCategory) -> List[str]:
        """Get all genes in a category"""
        genes = []
        for gene_name, info in cls.GENE_INFO.items():
            if info.category == category:
                genes.append(gene_name)
        return sorted(genes)
    
    @classmethod
    def get_drugs_by_gene(cls, gene: str) -> List[str]:
        """Get drugs with PD gene implications"""
        info = cls.get_gene_info(gene)
        if info:
            return info.drug_list
        return []
    
    @classmethod
    def search_genes(cls, search_term: str) -> List[PDGeneInfo]:
        """Search genes by name or function"""
        results = []
        search_lower = search_term.lower()
        
        for gene, info in cls.GENE_INFO.items():
            if (search_lower in gene.lower() or
                search_lower in info.full_name.lower() or
                search_lower in info.biological_function.lower() or
                search_lower in info.clinical_relevance.lower()):
                results.append(info)
        
        return results


class PDGeneRecommendationDatabase:
    """Recommendations for pharmacodynamic genes"""
    
    RECOMMENDATIONS = {
        ('HLA-B', '*5701', 'abacavir'): PDGeneRecommendation(
            gene='HLA-B',
            variant_or_phenotype='HLA-B*5701 positive',
            drug='abacavir',
            recommendation='AVOID - CONTRAINDICATED',
            clinical_rationale='5-8% hypersensitivity reaction risk; potentially fatal reaction',
            monitoring_parameters=['No use - contraindicated'],
            alternatives=['Integrase inhibitors (DTG, BIC)', 'Protease inhibitors', 'NNRTI combinations']
        ),
        
        ('HLA-A', '*31:01', 'carbamazepine'): PDGeneRecommendation(
            gene='HLA-A',
            variant_or_phenotype='HLA-A*31:01 positive',
            drug='carbamazepine',
            recommendation='Alternative or close monitoring',
            clinical_rationale='Increased SJS/TEN risk, particularly in European and Asian populations',
            monitoring_parameters=['Rash assessment at baseline and each visit', 'Fever monitoring', 'Mucosal symptoms'],
            alternatives=['Lamotrigine (genetic testing for HLA-DRB1 needed)', 'Levetiracetam', 'Valproate']
        ),
        
        ('F5', 'Leiden', 'warfarin'): PDGeneRecommendation(
            gene='F5',
            variant_or_phenotype='Factor V Leiden',
            drug='warfarin',
            recommendation='Use with caution; may need reduced dose',
            clinical_rationale='Leiden carriers have increased baseline thrombosis risk; combined with oral anticoagulation needs close monitoring',
            monitoring_parameters=['INR target 2-3', 'More frequent INR monitoring', 'Thromboembolism assessment'],
            alternatives=['DOAC (Dabigatran, Apixaban, Rivaroxaban)', 'Fondaparinux']
        ),
        
        ('F2', 'G20210A', 'warfarin'): PDGeneRecommendation(
            gene='F2',
            variant_or_phenotype='Prothrombin 20210A carrier',
            drug='warfarin',
            recommendation='Use standard dosing; pharmacogenomic prediction helpful',
            clinical_rationale='Increased baseline thrombosis risk; interact with warfarin metabolism',
            monitoring_parameters=['INR target 2-3', 'Weekly INR checks x4, then less frequent', 'Baseline thromboembolism risk'],
            alternatives=['DOAC if consistent INR difficult to maintain']
        ),
        
        ('MTHFR', 'T677T', 'methotrexate'): PDGeneRecommendation(
            gene='MTHFR',
            variant_or_phenotype='TT homozygote (35% enzyme activity)',
            drug='methotrexate',
            recommendation='May use with caution; increased toxicity monitoring',
            clinical_rationale='Reduced methylenetetrahydrofolate reductase activity may increase MTX toxicity',
            monitoring_parameters=['Increased folic acid supplementation (5mg daily)', 'CBC, LFTs weekly x4, then monthly', 'GI side effect monitoring'],
            alternatives=['Biologics (TNF inhibitors)', 'Other DMARDs']
        ),
        
        ('MC4R', 'Risk allele', 'olanzapine'): PDGeneRecommendation(
            gene='MC4R',
            variant_or_phenotype='MC4R obesity risk allele carrier',
            drug='olanzapine',
            recommendation='CAUTION: Increased weight gain risk',
            clinical_rationale='rs17782313 risk allele carriers have 2-3x higher risk of antipsychotic-induced weight gain',
            monitoring_parameters=['Baseline and monthly weight', 'Waist circumference', 'Metabolic panel (glucose, lipids)', 'Consider lifestyle intervention'],
            alternatives=['Aripiprazole', 'Ziprasidone', 'Lurasidone']
        ),
        
        ('TP53', 'Germline mutation', 'doxorubicin'): PDGeneRecommendation(
            gene='TP53',
            variant_or_phenotype='Li-Fraumeni syndrome (TP53 germline mutation)',
            drug='doxorubicin',
            recommendation='CAUTION: Avoid if possible; increased secondary cancer risk',
            clinical_rationale='TP53 mutation carriers at extremely high risk of treatment-related secondary malignancies',
            monitoring_parameters=['Cumulative dose limits (≤250 mg/m²)', 'Cardiotoxicity monitoring', 'Long-term cancer surveillance'],
            alternatives=['Taxanes', 'Platinum agents', 'Immunotherapy']
        ),
    }
    
    @classmethod
    def get_recommendation(cls, gene: str, variant_or_phenotype: str, drug: str) -> Optional[PDGeneRecommendation]:
        """Get recommendation for PD gene variant-drug pair"""
        key = (gene.upper(), variant_or_phenotype, drug.lower())
        return cls.RECOMMENDATIONS.get(key)
    
    @classmethod
    def format_recommendation(cls, recommendation: PDGeneRecommendation) -> str:
        """Format recommendation for clinical report"""
        lines = []
        lines.append(f"GENE: {recommendation.gene}")
        lines.append(f"VARIANT/PHENOTYPE: {recommendation.variant_or_phenotype}")
        lines.append(f"DRUG: {recommendation.drug}")
        lines.append(f"RECOMMENDATION: {recommendation.recommendation}")
        lines.append(f"\nRATIONALE: {recommendation.clinical_rationale}")
        
        if recommendation.monitoring_parameters:
            lines.append("\nMONITORING PARAMETERS:")
            for param in recommendation.monitoring_parameters:
                lines.append(f"  • {param}")
        
        if recommendation.alternatives:
            lines.append("\nALTERNATIVE DRUGS:")
            for alt in recommendation.alternatives:
                lines.append(f"  • {alt}")
        
        return "\n".join(lines)


def get_pharmacodynamic_genes_summary(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Generate summary of relevant pharmacodynamic genes
    
    Args:
        genotypes: Dictionary of detected gene genotypes/phenotypes
        
    Returns:
        Summary with clinical implications
    """
    summary = {
        'immune_risk_genes': [],
        'hemostasis_genes': [],
        'metabolic_risk_genes': [],
        'psychiatric_response_genes': [],
        'other_pd_genes': []
    }
    
    for gene, phenotype in genotypes.items():
        info = PharmacodynamicGenes.get_gene_info(gene)
        if info:
            entry = {
                'gene': gene,
                'phenotype': phenotype,
                'category': info.category.value,
                'clinical_relevance': info.clinical_relevance,
                'drugs_affected': info.drug_list
            }
            
            if info.category == GeneCategory.HYPERSENSITIVITY:
                summary['immune_risk_genes'].append(entry)
            elif info.category == GeneCategory.HEMOSTASIS:
                summary['hemostasis_genes'].append(entry)
            elif info.category == GeneCategory.METABOLIC:
                summary['metabolic_risk_genes'].append(entry)
            elif info.category == GeneCategory.PSYCHIATRIC:
                summary['psychiatric_response_genes'].append(entry)
            else:
                summary['other_pd_genes'].append(entry)
    
    return summary
