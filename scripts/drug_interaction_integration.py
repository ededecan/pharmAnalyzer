#!/usr/bin/env python3
"""
Drug Interaction Integration Module
Provides utilities for extracting and formatting drug interactions in NutriGeneEngine style
Can be used standalone or integrated into the main analysis pipeline

ENHANCED VERSION with:
- Risk stratification (RED/ORANGE/YELLOW/GREEN)
- CPIC guidelines integration
- Intelligent DDI filtering
- ATC therapeutic area classification
"""

# Configure matplotlib backend for non-interactive plotting (must be before pyplot import)
import matplotlib
matplotlib.use('Agg')

from typing import Dict, List, Any, Optional
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import re

def is_valid_drug_name(name: str) -> bool:
    """Validate drug name to filter out invalid entries like single chars or pure numbers"""
    if not name or not isinstance(name, str):
        return False
    name = name.strip()
    if len(name) <= 2:  # Filter 1-2 character names
        return False
    if name.isdigit():  # Filter pure numbers
        return False
    if not any(c.isalpha() for c in name):  # Must have at least one letter
        return False
    return True

# Import enhanced pharmacogenomics modules
try:
    from risk_stratification import RiskStratifier, RiskLevel, prioritize_interactions
    from cpic_guidelines import CPICGuidelines, CPICRecommendation
    from atc_classification import ATCClassifier, TherapeuticArea
    from ddi_filtering import DDIFilter, DDISeverity, summarize_ddi_statistics
    ENHANCED_FEATURES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Enhanced features not available: {e}")
    print("Continuing with basic functionality...")
    ENHANCED_FEATURES_AVAILABLE = False


class DrugBankXMLParser:
    """Streaming parser for DrugBank XML database"""
    
    @staticmethod
    def parse_drugbank_xml(xml_path: str) -> Dict[str, Any]:
        """Stream-parse DrugBank XML to extract interactions and targets
        
        Args:
            xml_path: Path to full-drug-database.xml
            
        Returns:
            Dict with 'interactions' and 'targets' lists
        """
        interactions = []
        targets = []
        
        try:
            # Use iterparse for memory-efficient streaming
            context = ET.iterparse(xml_path, events=('start', 'end'))
            context = iter(context)
            event, root = next(context)
            
            namespace = {'db': 'http://www.drugbank.ca'}
            current_drug = {}
            in_drug = False
            
            for event, elem in context:
                tag = elem.tag.replace('{http://www.drugbank.ca}', '')
                
                if event == 'start' and tag == 'drug':
                    in_drug = True
                    current_drug = {'name': '', 'drugbank_id': '', 'interactions': [], 'targets': []}
                    
                elif event == 'end' and in_drug:
                    if tag == 'drugbank-id' and elem.attrib.get('primary') == 'true':
                        current_drug['drugbank_id'] = elem.text or ''
                        
                    elif tag == 'name' and not current_drug['name']:
                        current_drug['name'] = elem.text or ''
                        
                    elif tag == 'drug-interaction':
                        interaction_data = {}
                        for child in elem:
                            child_tag = child.tag.replace('{http://www.drugbank.ca}', '')
                            if child_tag == 'drugbank-id':
                                interaction_data['interacting_drug_id'] = child.text
                            elif child_tag == 'name':
                                interaction_data['interacting_drug_name'] = child.text
                            elif child_tag == 'description':
                                interaction_data['description'] = child.text
                        if interaction_data:
                            current_drug['interactions'].append(interaction_data)
                            
                    elif tag == 'target':
                        target_data = {}
                        for child in elem:
                            child_tag = child.tag.replace('{http://www.drugbank.ca}', '')
                            if child_tag == 'polypeptide':
                                for poly_child in child:
                                    poly_tag = poly_child.tag.replace('{http://www.drugbank.ca}', '')
                                    if poly_tag == 'gene-name':
                                        target_data['gene'] = poly_child.text
                                    elif poly_tag == 'external-identifiers':
                                        for ext_id in poly_child:
                                            res_tag = ext_id.find('.//{http://www.drugbank.ca}resource')
                                            id_tag = ext_id.find('.//{http://www.drugbank.ca}identifier')
                                            if res_tag is not None and res_tag.text == 'UniProtKB':
                                                if id_tag is not None:
                                                    target_data['uniprot'] = id_tag.text
                        if target_data.get('gene'):
                            target_data['drug_name'] = current_drug['name']
                            target_data['drug_id'] = current_drug['drugbank_id']
                            current_drug['targets'].append(target_data)
                            
                    elif tag == 'drug':
                        # End of drug element - store data
                        for interaction in current_drug['interactions']:
                            interactions.append({
                                'query_drug': current_drug['name'],
                                'query_drug_id': current_drug['drugbank_id'],
                                'interacting_drug_id': interaction.get('interacting_drug_id', ''),
                                'interacting_drug_name': interaction.get('interacting_drug_name', ''),
                                'description': interaction.get('description', ''),
                                'severity': DrugBankXMLParser._infer_severity(interaction.get('description', ''))
                            })
                        
                        for target in current_drug['targets']:
                            targets.append(target)
                        
                        in_drug = False
                        current_drug = {}
                        elem.clear()
                        root.clear()
                        
        except Exception as e:
            print(f"Warning: XML parsing error: {e}")
            import traceback
            traceback.print_exc()
            
        return {'interactions': interactions, 'targets': targets}
    
    @staticmethod
    def _infer_severity(description: str) -> str:
        """Infer interaction severity from description text"""
        if ENHANCED_FEATURES_AVAILABLE:
            severity = DDIFilter.infer_severity(description)
            return severity.value
        return 'UNKNOWN'


class DrugInteractionExtractor:
    """Extracts drug interaction data from analysis results"""
    
    @staticmethod
    def extract_from_analysis(analysis_results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract all drug interactions from analysis results"""
        interactions = {
            'sample_id': analysis_results.get('sample_id', 'Unknown'),
            'timestamp': analysis_results.get('analysis_date', ''),
            'genes_analyzed': {},
            'critical_interactions': [],
            'moderate_interactions': [],
            'general_interactions': []
        }
        
        for gene, gene_data in analysis_results.get('variants_by_gene', {}).items():
            diplotype = gene_data.get('diplotype', 'Unknown')
            
            gene_info = {
                'gene': gene,
                'diplotype': diplotype,
                'phenotype': gene_data.get('coded_summary', 'Unknown'),
                'metabolizer_status': gene_data.get('metabolizer_status', 'Unknown'),
                'drugs': {}
            }
            
            for rec in gene_data.get('pharmgkb_recommendations', []):
                drug = rec.get('drug', 'Unknown')
                action = rec.get('action', 'CONSULT SPECIALIST').upper()
                phenotype = gene_data.get('metabolizer_status', 'Unknown')
                evidence_level = rec.get('evidence_level', '')
                
                interaction_data = {
                    'gene': gene,
                    'diplotype': diplotype,
                    'drug': drug,
                    'action': action,
                    'direction': rec.get('direction', ''),
                    'reason': rec.get('reason', ''),
                    'pmid': rec.get('pmid', ''),
                    'notes': rec.get('notes', ''),
                    'evidence_level': evidence_level
                }
                
                if ENHANCED_FEATURES_AVAILABLE:
                    risk_level = RiskStratifier.calculate_risk(
                        gene, diplotype, phenotype, drug, evidence_level
                    )
                    interaction_data['risk_level'] = risk_level.value
                    interaction_data['risk_color'] = risk_level.get_color_hex()
                    interaction_data['risk_description'] = risk_level.get_description()
                    
                    action_rec = RiskStratifier.get_action_recommendation(
                        risk_level, gene, phenotype, drug
                    )
                    interaction_data['action_recommendation'] = action_rec
                    
                    warning = RiskStratifier.get_critical_warning(
                        gene, diplotype, phenotype, drug
                    )
                    if warning:
                        interaction_data['critical_warning'] = warning
                    
                    cpic_rec = CPICGuidelines.get_recommendation(gene, phenotype, drug)
                    if cpic_rec:
                        interaction_data['cpic_guideline'] = cpic_rec.to_dict()
                    
                    therapeutic_area = ATCClassifier.get_therapeutic_area(drug)
                    interaction_data['therapeutic_area'] = therapeutic_area.value
                
                if action == 'AVOID':
                    interactions['critical_interactions'].append(interaction_data)
                elif action == 'ADJUST DOSE':
                    interactions['moderate_interactions'].append(interaction_data)
                else:
                    interactions['general_interactions'].append(interaction_data)
                
                gene_info['drugs'][drug] = {
                    'action': action,
                    'direction': rec.get('direction', ''),
                    'notes': rec.get('notes', '')
                }
            
            interactions['genes_analyzed'][gene] = gene_info
        
        interactions['summary'] = {
            'total_genes': len(interactions['genes_analyzed']),
            'critical_count': len(interactions['critical_interactions']),
            'moderate_count': len(interactions['moderate_interactions']),
            'general_count': len(interactions['general_interactions']),
            'total_interactions': (len(interactions['critical_interactions']) + 
                                  len(interactions['moderate_interactions']) + 
                                  len(interactions['general_interactions']))
        }
        
        return interactions
    
    @staticmethod
    def get_critical_interactions(interactions: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Get all critical (AVOID) interactions
        
        Args:
            interactions: Extracted interaction data
        
        Returns:
            List of critical interactions
        """
        return interactions.get('critical_interactions', [])
    
    @staticmethod
    def get_interactions_for_drug(interactions: Dict[str, Any], drug: str) -> List[Dict[str, Any]]:
        """Get all interactions for a specific drug
        
        Args:
            interactions: Extracted interaction data
            drug: Drug name (case-insensitive)
        
        Returns:
            List of interactions for this drug
        """
        drug_lower = drug.lower()
        all_interactions = (interactions.get('critical_interactions', []) + 
                           interactions.get('moderate_interactions', []) + 
                           interactions.get('general_interactions', []))
        
        return [i for i in all_interactions if i.get('drug', '').lower() == drug_lower]
    
    @staticmethod
    def get_interactions_for_gene(interactions: Dict[str, Any], gene: str) -> List[Dict[str, Any]]:
        """Get all interactions for a specific gene
        
        Args:
            interactions: Extracted interaction data
            gene: Gene name (case-insensitive)
        
        Returns:
            List of interactions for this gene
        """
        gene_lower = gene.lower()
        all_interactions = (interactions.get('critical_interactions', []) + 
                           interactions.get('moderate_interactions', []) + 
                           interactions.get('general_interactions', []))
        
        return [i for i in all_interactions if i.get('gene', '').lower() == gene_lower]


def integrate_drug_interactions_into_main_analysis(results: Dict[str, Any],
                                                   patient_meds: List[str],
                                                   db_content: Optional[str] = None) -> Dict[str, Any]:
    """Integrate drug interaction information into the main analysis results.

    This function extracts interaction-level recommendations from the analysis
    results (for example, PharmGKB-derived recommendations), and optionally
    parses DrugBank XML to extract drug-drug interactions and drug-target mappings.
    It then matches recommendations and interactions against the patient's current
    medication list and attaches a `clinical_decision_support` block to the results dict.

    Args:
        results: The analysis results dictionary produced by the analyzer.
        patient_meds: List of patient medication names (may be empty).
        db_content: Optional path to a drug DB (XML path) or text content.

    Returns:
        The updated `results` dict with a `clinical_decision_support` key added.
    """
    print("\n🔬 Integrating drug interaction data...")
    
    # Extract structured interactions from PharmGKB analysis results
    interactions = DrugInteractionExtractor.extract_from_analysis(results)
    
    # Prepare patient meds lowercase set for matching
    meds_set = set(m.strip().lower() for m in (patient_meds or []) if m and m.strip())
    
    # Parse DrugBank XML if provided
    drugbank_data = {'interactions': [], 'targets': []}
    db_source = 'PharmGKB/Local'
    
    if db_content and isinstance(db_content, str):
        if db_content.lower().endswith('.xml') and Path(db_content).exists():
            print(f"📊 Parsing DrugBank XML: {Path(db_content).name}")
            print("   This may take several minutes for large files...")
            drugbank_data = DrugBankXMLParser.parse_drugbank_xml(db_content)
            print(f"   ✓ Extracted {len(drugbank_data['interactions'])} drug interactions")
            print(f"   ✓ Extracted {len(drugbank_data['targets'])} drug-target mappings")
            db_source = 'DrugBank (XML) + PharmGKB'
    
    gene_to_drugbank_drugs = {}
    for target in drugbank_data['targets']:
        gene = target.get('gene', '').upper()
        if gene and gene in results.get('variants_by_gene', {}):
            if gene not in gene_to_drugbank_drugs:
                gene_to_drugbank_drugs[gene] = []
            gene_to_drugbank_drugs[gene].append({
                'drug_name': target.get('drug_name', ''),
                'drug_id': target.get('drug_id', ''),
                'uniprot': target.get('uniprot', ''),
                'source': 'DrugBank'
            })
    
    print(f"   📋 Building gene-specific drug lists from PharmGKB...")
    gene_pharmgkb_drugs = {}
    for gene, gene_data in results.get('variants_by_gene', {}).items():
        drugs_for_gene = set()
        for drug in gene_data.get('affected_drugs', []):
            if drug and drug.strip() and is_valid_drug_name(drug.strip()):
                drugs_for_gene.add(drug.strip().lower())
        for rec in gene_data.get('pharmgkb_recommendations', []):
            drug = rec.get('drug', '')
            if drug and drug.strip() and is_valid_drug_name(drug.strip()):
                drugs_for_gene.add(drug.strip().lower())
        
        if drugs_for_gene:
            gene_pharmgkb_drugs[gene] = drugs_for_gene
    
    gene_drugbank_ddis = {}
    if drugbank_data['interactions'] and gene_pharmgkb_drugs:
        print(f"   🔍 Filtering DrugBank DDIs for gene-associated drugs...")
        
        for gene, gene_drugs in gene_pharmgkb_drugs.items():
            # Track unique interactions per drug pair: (drug1, drug2) -> list of interactions
            gene_ddis_by_pair = {}  
            
            for ddi in drugbank_data['interactions']:
                query_drug = (ddi.get('query_drug') or '').strip().lower()
                interacting_drug = (ddi.get('interacting_drug_name') or '').strip().lower()
                description = (ddi.get('description') or '').strip()
                
                # Include DDI if either drug is associated with this gene in PharmGKB
                if query_drug in gene_drugs or interacting_drug in gene_drugs:
                    # Use ordered pair (don't sort) to preserve directionality
                    drug_pair = (query_drug, interacting_drug)
                    
                    if drug_pair not in gene_ddis_by_pair:
                        gene_ddis_by_pair[drug_pair] = []
                    
                    # Add this interaction (deduplicate by description to avoid exact duplicates)
                    interaction_entry = {
                        'drug_1': ddi.get('query_drug', ''),
                        'drug_2': ddi.get('interacting_drug_name', ''),
                        'description': description,
                        'severity': ddi.get('severity', 'UNKNOWN'),
                        'source': 'DrugBank'
                    }
                    
                    # Only add if this exact description doesn't already exist for this pair
                    desc_exists = any(
                        entry.get('description', '').strip().lower() == description.lower() 
                        for entry in gene_ddis_by_pair[drug_pair]
                    )
                    if not desc_exists:
                        gene_ddis_by_pair[drug_pair].append(interaction_entry)
            
            # Flatten to list - keeping all unique interaction types
            if gene_ddis_by_pair:
                gene_ddis_list = []
                for pair_interactions in gene_ddis_by_pair.values():
                    gene_ddis_list.extend(pair_interactions)
                
                # Apply intelligent filtering if enhanced features are available
                if ENHANCED_FEATURES_AVAILABLE:
                    print(f"      🧹 Applying intelligent filtering for {gene} ({len(gene_ddis_list)} raw DDIs)...")
                    
                    # Keep a copy of the original unfiltered list for fallback passes
                    original_unfiltered = list(gene_ddis_list)
                    filtered_result = []
                    
                    # Pass 1: Strict filtering (MODERATE+ severity, common drugs, gene-specific)
                    filtered_result = DDIFilter.filter_ddis(
                        original_unfiltered,
                        gene=gene,
                        severity_threshold=DDISeverity.MODERATE,
                        common_drugs_only=True,
                        gene_specific_only=True,
                        max_interactions=200
                    )
                    
                    # Pass 2: If too few results, relax gene-specific requirement
                    if len(filtered_result) < 30:
                        filtered_result = DDIFilter.filter_ddis(
                            original_unfiltered,
                            gene=gene,
                            severity_threshold=DDISeverity.MODERATE,
                            common_drugs_only=True,
                            gene_specific_only=False,  # Allow all drug interactions
                            max_interactions=200
                        )
                    
                    # Pass 3: If still too few, relax common drugs requirement
                    if len(filtered_result) < 30:
                        filtered_result = DDIFilter.filter_ddis(
                            original_unfiltered,
                            gene=gene,
                            severity_threshold=DDISeverity.MODERATE,
                            common_drugs_only=False,  # Allow any drug
                            gene_specific_only=False,
                            max_interactions=200
                        )
                    
                    # Pass 4: If still too few, lower severity threshold
                    if len(filtered_result) < 30:
                        filtered_result = DDIFilter.filter_ddis(
                            original_unfiltered,
                            gene=gene,
                            severity_threshold=DDISeverity.MINOR,  # Include MINOR severity
                            common_drugs_only=False,
                            gene_specific_only=False,
                            max_interactions=200
                        )
                    
                    # Deduplicate
                    gene_ddis_list = DDIFilter.deduplicate_ddis(filtered_result)
                    
                    # Ensure we have results but cap at 200 for performance
                    if len(gene_ddis_list) > 200:
                        gene_ddis_list = gene_ddis_list[:200]
                    
                    print(f"      ✓ Filtered to {len(gene_ddis_list)} clinically relevant DDIs")
                
                gene_drugbank_ddis[gene] = gene_ddis_list
        
        total_gene_ddis = sum(len(ddis) for ddis in gene_drugbank_ddis.values())
        print(f"   ✓ Found {total_gene_ddis} unique gene-specific DrugBank DDI entries across {len(gene_drugbank_ddis)} genes")
        if ENHANCED_FEATURES_AVAILABLE:
            print(f"      (Filtered for clinical relevance: MODERATE+ severity, common drugs, gene-specific substrates)")
        else:
            print(f"      (Counting unique interaction types per drug pair for frequency-based risk assessment)")
    
    # Find PharmGKB-based alerts where the recommended drug is in patient's meds
    interaction_alerts: List[Dict[str, Any]] = []
    
    for category in ('critical_interactions', 'moderate_interactions', 'general_interactions'):
        for intr in interactions.get(category, []):
            drug_name = (intr.get('drug') or '').strip().lower()
            if drug_name and drug_name in meds_set:
                alert = dict(intr)
                alert['source'] = 'PharmGKB'
                alert['_matched_category'] = category
                interaction_alerts.append(alert)
    
    # Add patient medication-based DrugBank DDI alerts (for clinical decision support view)
    if meds_set and drugbank_data['interactions']:
        print(f"   🔍 Checking DrugBank interactions for {len(meds_set)} patient medications...")
        ddi_count = 0
        
        for ddi in drugbank_data['interactions']:
            query_drug = (ddi.get('query_drug') or '').strip().lower()
            interacting_drug = (ddi.get('interacting_drug_name') or '').strip().lower()
            
            # Check if either drug is in patient's meds
            if query_drug in meds_set or interacting_drug in meds_set:
                interaction_alerts.append({
                    'drug_1': ddi.get('query_drug', ''),
                    'drug_2': ddi.get('interacting_drug_name', ''),
                    'description': ddi.get('description', ''),
                    'severity': ddi.get('severity', 'UNKNOWN'),
                    'source': 'DrugBank',
                    '_matched_category': 'drugbank_ddi'
                })
                ddi_count += 1
        
        if ddi_count > 0:
            print(f"   ✓ Found {ddi_count} DrugBank drug-drug interactions")
    
    # Merge gene-drug annotations from PharmGKB and DrugBank
    gene_drug_annotations = interactions.get('genes_analyzed', {})
    
    for gene, gene_data in gene_drug_annotations.items():
        # Add DrugBank targets if available
        if gene in gene_to_drugbank_drugs:
            gene_data['drugbank_targets'] = gene_to_drugbank_drugs[gene]
        # Add gene-specific DrugBank DDIs
        if gene in gene_drugbank_ddis:
            gene_data['drugbank_ddis'] = gene_drugbank_ddis[gene]
    
    # Also attach gene-specific DrugBank DDIs directly to variants_by_gene for PDF rendering
    for gene in gene_drugbank_ddis:
        if gene in results.get('variants_by_gene', {}):
            results['variants_by_gene'][gene]['drugbank_ddis'] = gene_drugbank_ddis[gene]
    
    # Build the clinical decision support block
    cds_block: Dict[str, Any] = {
        'patient_medications': sorted(list(meds_set)),
        'interaction_alerts': interaction_alerts,
        'gene_drug_annotations': gene_drug_annotations,
        'interaction_database_source': db_source,
        'drugbank_summary': {
            'total_interactions': len(drugbank_data['interactions']),
            'total_targets': len(drugbank_data['targets']),
            'genes_with_drugbank_targets': len(gene_to_drugbank_drugs),
            'genes_with_drugbank_ddis': len(gene_drugbank_ddis)
        }
    }
    
    results['clinical_decision_support'] = cds_block
    
    print(f"✓ Clinical Decision Support integrated:")
    print(f"   • {len(meds_set)} patient medications")
    print(f"   • {len(interaction_alerts)} interaction alerts")
    print(f"   • {len(gene_drug_annotations)} genes annotated")
    if gene_drugbank_ddis:
        print(f"   • {len(gene_drugbank_ddis)} genes with DrugBank DDIs")
    
    return results


class DDIHeatmapGenerator:
    """Generate heatmap visualizations of drug-drug interactions"""
    
    @staticmethod
    def generate_gene_interaction_heatmap(gene: str, ddis: List[Dict[str, Any]], 
                                         affected_drugs: List[str],
                                         output_path: str) -> None:
        """Generate a heatmap for a specific gene showing ALL drug interactions
        
        Approach: Use DrugBank severity for each drug pair (worst-case severity)
        - Y axis: Drugs affected by gene (from PharmGKB)
        - X axis: All drugs that interact with Y-axis drugs (from DrugBank)
        - Color intensity: Severity level (LOW/MODERATE/HIGH)
        
        Args:
            gene: Gene name
            ddis: List of DDI dictionaries from DrugBank XML
            affected_drugs: List of drugs affected by this gene from PharmGKB
            output_path: Path to save the heatmap PNG file
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            import pandas as pd
            from matplotlib.colors import ListedColormap
        except ImportError as e:
            print(f"⚠️  Warning: Required libraries not available: {e}")
            return
        
        if not ddis or not affected_drugs:
            return
        
        # Y-axis: Split compound drug names from PharmGKB and clean (rows)
        affected_drugs_clean = []
        for drug_entry in affected_drugs:
            if not drug_entry or not drug_entry.strip():
                continue
            # Split by comma to handle "drug1, drug2" entries
            drugs_list = [d.strip() for d in drug_entry.split(',')]
            for drug in drugs_list:
                # Clean up and standardize
                drug_clean = drug.strip().title()  # Capitalize for matching with DrugBank
                if drug_clean and drug_clean not in affected_drugs_clean and is_valid_drug_name(drug_clean):
                    affected_drugs_clean.append(drug_clean)
        
        if not affected_drugs_clean:
            print(f"⚠ Warning: No affected drugs for {gene}")
            return
        
        # Build interaction severity map: (drug1, drug2) -> max severity score
        severity_rank = {
            'SEVERE': 3,
            'MAJOR': 3,
            'HIGH': 3,
            'MODERATE': 2,
            'MEDIUM': 2,
            'MINOR': 1,
            'LOW': 1,
            'UNKNOWN': 0
        }
        interaction_severity = {}
        
        for ddi in ddis:
            drug1 = ddi.get('drug_1', '').strip().title()
            drug2 = ddi.get('drug_2', '').strip().title()
            
            if not drug1 or not drug2 or not is_valid_drug_name(drug1) or not is_valid_drug_name(drug2):
                continue
            
            # Check if either drug is in our affected drugs list
            drug1_affected = any(drug1.lower() == affected.lower() for affected in affected_drugs_clean)
            drug2_affected = any(drug2.lower() == affected.lower() for affected in affected_drugs_clean)
            
            if not (drug1_affected or drug2_affected):
                continue
            
            severity_raw = (ddi.get('severity') or 'UNKNOWN').strip().upper()
            severity_score = severity_rank.get(severity_raw, 0)
            
            # Store as (affected_drug, interacting_drug) for directed interaction
            if drug1_affected:
                key = (drug1, drug2)
                interaction_severity[key] = max(interaction_severity.get(key, 0), severity_score)
            if drug2_affected and drug1 != drug2:
                key = (drug2, drug1)
                interaction_severity[key] = max(interaction_severity.get(key, 0), severity_score)
        
        if not interaction_severity:
            print(f"⚠ Warning: No matching drugs found for {gene} in DDI database")
            return
        
        # X-axis: Get all unique interacting drugs (drugs that interact with affected drugs)
        all_interacting_drugs = set()
        for (drug1, drug2) in interaction_severity.keys():
            if any(drug1.lower() == affected.lower() for affected in affected_drugs_clean):
                all_interacting_drugs.add(drug2)
            if any(drug2.lower() == affected.lower() for affected in affected_drugs_clean):
                all_interacting_drugs.add(drug1)
        
        # Count total interactions per X-axis drug to sort by most interactions
        x_drug_totals = {}
        for x_drug in all_interacting_drugs:
            count = 0
            for (drug1, drug2), score in interaction_severity.items():
                if x_drug.lower() == drug2.lower() or x_drug.lower() == drug1.lower():
                    count += score
            x_drug_totals[x_drug] = count
        
        # Sort X-axis drugs by total interaction count (descending)
        sorted_x_drugs = sorted(x_drug_totals.items(), key=lambda x: x[1], reverse=True)
        
        # Limit to top 60 X-axis drugs for readability
        if len(sorted_x_drugs) > 60:
            print(f"⚠ Note: Limiting {gene} heatmap to top 60 interacting drugs (from {len(sorted_x_drugs)})")
            sorted_x_drugs = sorted_x_drugs[:60]
        
        all_interacting_drugs = [d[0] for d in sorted_x_drugs]
        
        # Y-axis: Sort affected drugs by total interaction count and limit to top 60 clinically significant drugs
        y_drug_totals = {}
        for y_drug in affected_drugs_clean:
            count = 0
            for (drug1, drug2), score in interaction_severity.items():
                if y_drug.lower() == drug1.lower() or y_drug.lower() == drug2.lower():
                    count += score
            y_drug_totals[y_drug] = count
        
        # Sort by total interactions (descending) to show most clinically significant drugs
        sorted_y_drugs = sorted(y_drug_totals.items(), key=lambda x: x[1], reverse=True)
        
        # Limit to top 60 Y-axis drugs for clinical relevance
        if len(sorted_y_drugs) > 60:
            print(f"⚠ Note: Limiting {gene} heatmap to top 60 clinically significant affected drugs (from {len(sorted_y_drugs)})")
            sorted_y_drugs = sorted_y_drugs[:60]
        
        affected_drugs_clean = [d[0] for d in sorted_y_drugs]
        
        if len(affected_drugs_clean) < 1 or len(all_interacting_drugs) < 1:
            return
        
        # Create matrix
        n_rows = len(affected_drugs_clean)
        n_cols = len(all_interacting_drugs)
        matrix = np.zeros((n_rows, n_cols))
        
        # Fill matrix with interaction severity (worst case per pair)
        for i, y_drug in enumerate(affected_drugs_clean):
            for j, x_drug in enumerate(all_interacting_drugs):
                score = 0
                for (drug1, drug2), severity_score in interaction_severity.items():
                    if ((y_drug.lower() == drug1.lower() and x_drug.lower() == drug2.lower()) or
                        (y_drug.lower() == drug2.lower() and x_drug.lower() == drug1.lower())):
                        score = max(score, severity_score)
                matrix[i, j] = score
        
        # Create DataFrame
        df = pd.DataFrame(matrix, index=affected_drugs_clean, columns=all_interacting_drugs)
        
        # Get matrix statistics for color scaling
        max_score = int(matrix.max())
        min_score = int(matrix[matrix > 0].min()) if np.any(matrix > 0) else 0
        
        print(f"📊 {gene} interaction matrix: {n_rows}×{n_cols}, "
              f"max severity: {max_score}, min: {min_score}")
        
        # Create figure
        try:
            # Calculate dimensions based on filtered data
            n_rows = len(affected_drugs_clean)
            n_cols = len(all_interacting_drugs)
            
            # Set reasonable figure size
            figsize = (min(24, max(14, n_cols * 0.4)), min(16, max(8, n_rows * 0.5)))
            fig, ax = plt.subplots(figsize=figsize, dpi=100)
            
            # Define color scheme based on SEVERITY of interactions
            from matplotlib.colors import ListedColormap, BoundaryNorm
            import matplotlib.patches as mpatches
            
            # Map severity scores to severity categories
            colors = ['#f0f0f0', '#FFD700', '#FF8C00', '#DC143C']  # gray, yellow, orange, red
            boundaries = [0, 0.5, 1.5, 2.5, 3.5]
            color_list = colors
            labels = ['None/No Data', 'LOW', 'MODERATE', 'HIGH']
            
            # Create discrete colormap with proper normalization
            cmap = ListedColormap(color_list)
            norm = BoundaryNorm(boundaries, cmap.N)
            
            # Generate heatmap with discrete colors
            im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect='auto')
            
            # Add text annotations showing severity levels
            for i in range(n_rows):
                for j in range(n_cols):
                    score = int(matrix[i, j])
                    if score <= 0:
                        continue
                    if score == 3:
                        ax.text(j, i, 'H',
                               ha='center', va='center',
                               fontsize=10, fontweight='bold', color='white')
                    elif score == 2:
                        ax.text(j, i, 'M',
                               ha='center', va='center',
                               fontsize=9, fontweight='normal', color='white')
                    elif score == 1:
                        ax.text(j, i, 'L',
                               ha='center', va='center',
                               fontsize=8, fontweight='normal', color='black')
            
            # Shorten drug names for readability (clean truncation without ...)
            short_row_labels = [d[:25] if len(d) > 25 else d for d in affected_drugs_clean]
            short_col_labels = [d[:25] if len(d) > 25 else d for d in all_interacting_drugs]
            
            # Set axis labels
            ax.set_xticks(np.arange(n_cols))
            ax.set_yticks(np.arange(n_rows))
            ax.set_xticklabels(short_col_labels, rotation=90, ha='right', fontsize=9)
            ax.set_yticklabels(short_row_labels, fontsize=10, rotation=0, ha='right')
            
            ax.set_xlabel('Drug 2 (Interacting)', fontsize=13, fontweight='bold', labelpad=10)
            ax.set_ylabel('Drug 1 (Gene-Affected)', fontsize=13, fontweight='bold', labelpad=10)
            
            # Title with clear explanation
            title_text = (f'{gene} Drug-Drug Interaction Matrix\n'
                         f'{n_rows} gene-affected drugs × {n_cols} interacting drugs')
            ax.set_title(title_text, fontsize=14, fontweight='bold', pad=28)

            # Severity legend on right side (vertical)
            legend_elements = [
                mpatches.Patch(facecolor='#DC143C', edgecolor='#DC143C', label='HIGH (H)'),
                mpatches.Patch(facecolor='#FF8C00', edgecolor='#FF8C00', label='MODERATE (M)'),
                mpatches.Patch(facecolor='#2E8B57', edgecolor='#2E8B57', label='LOW (L)')
            ]
            legend = ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5),
                             fontsize=12, frameon=True, fancybox=False, edgecolor='black',
                             title='Severity', title_fontsize=13, framealpha=0.95)
            
            # Add white grid for cell separation
            ax.set_xticks(np.arange(n_cols + 1) - 0.5, minor=True)
            ax.set_yticks(np.arange(n_rows + 1) - 0.5, minor=True)
            ax.grid(which='minor', color='white', linestyle='-', linewidth=2)
            ax.tick_params(which='minor', size=0)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
            plt.close()
            
        except Exception as e:
            print(f"Error generating heatmap for {gene}: {e}")
            import traceback
            traceback.print_exc()
            return
    
    @staticmethod
    def generate_all_gene_heatmaps(results: Dict[str, Any], output_prefix: str) -> Dict[str, str]:
        """Generate interaction heatmaps for all genes
        
        Args:
            results: Analysis results with clinical_decision_support and variants_by_gene
            output_prefix: Prefix for output files
            
        Returns:
            Dict mapping gene names to heatmap file paths
        """
        cds = results.get('clinical_decision_support', {})
        gene_annotations = cds.get('gene_drug_annotations', {})
        variants_by_gene = results.get('variants_by_gene', {})
        
        gene_heatmap_paths = {}
        
        for gene, gene_data in gene_annotations.items():
            ddis = gene_data.get('drugbank_ddis', [])
            
            # Get affected drugs for this gene
            gene_variant_data = variants_by_gene.get(gene, {})
            affected_drugs = gene_variant_data.get('affected_drugs', [])
            
            if ddis:
                output_path = f"{output_prefix}_heatmap_{gene}.png"
                DDIHeatmapGenerator.generate_gene_interaction_heatmap(
                    gene, ddis, affected_drugs, output_path
                )
                gene_heatmap_paths[gene] = output_path
        
        if gene_heatmap_paths:
            print(f"✓ Generated {len(gene_heatmap_paths)} gene-specific interaction heatmaps")
        
        return gene_heatmap_paths
    
    @staticmethod
    def generate_interaction_heatmap(results: Dict[str, Any], output_path: str) -> None:
        """Generate a heatmap showing drug-drug interaction patterns across genes
        
        Args:
            results: Analysis results with clinical_decision_support
            output_path: Path to save the heatmap PNG file
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
            import pandas as pd
        except ImportError:
            print("⚠️  Warning: matplotlib/seaborn not installed. Skipping heatmap generation.")
            print("   Install with: pip install matplotlib seaborn")
            return
        
        cds = results.get('clinical_decision_support', {})
        gene_annotations = cds.get('gene_drug_annotations', {})
        
        # Collect all drugs and their interaction counts per gene
        drug_interaction_matrix = {}
        all_drugs = set()
        
        for gene, gene_data in gene_annotations.items():
            ddis = gene_data.get('drugbank_ddis', [])
            if not ddis:
                continue
            
            # Count interactions per drug for this gene
            drug_counts = {}
            for ddi in ddis:
                drug1 = ddi.get('drug_1', '').strip()
                drug2 = ddi.get('drug_2', '').strip()
                severity = ddi.get('severity', 'UNKNOWN')
                
                # Weight by severity
                weight = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'UNKNOWN': 0.5}.get(severity, 0.5)
                
                for drug in [drug1, drug2]:
                    if drug:
                        all_drugs.add(drug)
                        drug_counts[drug] = drug_counts.get(drug, 0) + weight
            
            drug_interaction_matrix[gene] = drug_counts
        
        if not drug_interaction_matrix:
            print("⚠️  No DDI data available for heatmap generation")
            return
        
        # Instead of gene × drug, create a drug × drug interaction matrix
        # This makes it clearer which drug pairs should not be combined
        
        # Collect all drug pairs and their highest severity
        all_drug_pairs = {}  # (drug1, drug2) -> max_severity_weight
        severity_weight = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'UNKNOWN': 0.5}
        
        for gene, gene_data in gene_annotations.items():
            ddis = gene_data.get('drugbank_ddis', [])
            for ddi in ddis:
                drug1 = ddi.get('drug_1', '').strip()
                drug2 = ddi.get('drug_2', '').strip()
                severity = ddi.get('severity', 'UNKNOWN')
                
                if drug1 and drug2:
                    key = tuple(sorted([drug1, drug2]))
                    weight = severity_weight.get(severity, 0.5)
                    
                    if key not in all_drug_pairs:
                        all_drug_pairs[key] = weight
                    else:
                        all_drug_pairs[key] = max(all_drug_pairs[key], weight)
        
        # Get top interacting drugs (by total interaction count)
        drug_totals = {}
        for (d1, d2), weight in all_drug_pairs.items():
            drug_totals[d1] = drug_totals.get(d1, 0) + weight
            drug_totals[d2] = drug_totals.get(d2, 0) + weight
        
        top_drugs = sorted(drug_totals.keys(), key=lambda d: drug_totals[d], reverse=True)[:25]
        
        # Build symmetric drug-drug matrix
        n = len(top_drugs)
        matrix_data = np.zeros((n, n))
        
        for i, drug1 in enumerate(top_drugs):
            for j, drug2 in enumerate(top_drugs):
                if i != j:
                    key = tuple(sorted([drug1, drug2]))
                    if key in all_drug_pairs:
                        matrix_data[i, j] = all_drug_pairs[key]
        
        # Create DataFrame
        df = pd.DataFrame(matrix_data, index=top_drugs, columns=top_drugs)
        
        # Generate improved heatmap
        from matplotlib.colors import LinearSegmentedColormap
        
        plt.figure(figsize=(18, 16))
        
        # Color scheme: light gray -> yellow -> orange -> red
        colors_list = [
            (0.0, '#f5f5f5'),   # Light gray for no interaction
            (0.33, '#ffeb3b'),  # Yellow for LOW
            (0.67, '#ff9800'),  # Orange for MODERATE
            (1.0, '#c62828')    # Dark red for HIGH
        ]
        cmap = LinearSegmentedColormap.from_list('severity', [c for _, c in colors_list], N=256)
        
        # Create heatmap with annotations for HIGH severity
        ax = sns.heatmap(df, cmap=cmap, annot=False, square=True,
                    cbar_kws={'label': 'Interaction Severity', 'shrink': 0.6,
                             'ticks': [0, 1, 2, 3]},
                    linewidths=0.8, linecolor='white', vmin=0, vmax=3)
        
        # Add warning symbols for HIGH severity pairs
        for i in range(n):
            for j in range(n):
                if matrix_data[i, j] == 3:
                    ax.text(j + 0.5, i + 0.5, '⚠', ha='center', va='center',
                           fontsize=12, fontweight='bold', color='white')
        
        # Customize colorbar
        cbar = ax.collections[0].colorbar
        cbar.set_ticklabels(['None', 'LOW', 'MODERATE', 'HIGH'])
        
        # Shorten drug names        
        short_labels = [d[:18] + '...' if len(d) > 18 else d for d in top_drugs]
        
        plt.title('Drug-Drug Interaction Matrix (All Genes Combined)\n'
                  '⚠ = HIGH Risk - DO NOT COMBINE THESE DRUGS\n'
                  f'Top {n} Most Frequently Interacting Drugs', 
                  fontsize=14, fontweight='bold', pad=20)
        plt.xlabel('Drug 2', fontsize=12, fontweight='bold')
        plt.ylabel('Drug 1', fontsize=12, fontweight='bold')
        plt.xticks(np.arange(n) + 0.5, short_labels, rotation=45, ha='right', fontsize=8)
        plt.yticks(np.arange(n) + 0.5, short_labels, rotation=0, fontsize=8)
        plt.tight_layout()
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ DDI heatmap saved to: {output_path}")
    
    @staticmethod
    def generate_severity_distribution(results: Dict[str, Any], output_path: str) -> None:
        """Generate a bar chart showing DDI severity distribution across genes
        
        Args:
            results: Analysis results with clinical_decision_support
            output_path: Path to save the chart PNG file
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            print("⚠️  Warning: matplotlib not installed. Skipping chart generation.")
            return
        
        cds = results.get('clinical_decision_support', {})
        gene_annotations = cds.get('gene_drug_annotations', {})
        
        # Collect severity counts per gene
        gene_severity_data = {}
        
        for gene, gene_data in gene_annotations.items():
            ddis = gene_data.get('drugbank_ddis', [])
            if not ddis:
                continue
            
            severity_counts = {'HIGH': 0, 'MODERATE': 0, 'LOW': 0}
            for ddi in ddis:
                severity = ddi.get('severity', 'UNKNOWN')
                if severity in severity_counts:
                    severity_counts[severity] += 1
            
            gene_severity_data[gene] = severity_counts
        
        if not gene_severity_data:
            print("⚠️  No DDI data available for severity chart")
            return
        
        # Prepare data for stacked bar chart
        genes = sorted(gene_severity_data.keys())
        high = [gene_severity_data[g]['HIGH'] for g in genes]
        moderate = [gene_severity_data[g]['MODERATE'] for g in genes]
        low = [gene_severity_data[g]['LOW'] for g in genes]
        
        # Create chart
        x = np.arange(len(genes))
        width = 0.6
        
        fig, ax = plt.subplots(figsize=(14, 7))
        
        # Stacked bars with clearer colors matching the heatmap
        p1 = ax.bar(x, high, width, label='⚠ HIGH (DO NOT COMBINE)', color='#c62828', edgecolor='white', linewidth=1.5)
        p2 = ax.bar(x, moderate, width, bottom=high, label='! MODERATE (Use Caution)', color='#ff9800', edgecolor='white', linewidth=1.5)
        p3 = ax.bar(x, low, width, bottom=np.array(high)+np.array(moderate), 
               label='LOW (Monitor)', color='#ffeb3b', edgecolor='white', linewidth=1.5)
        
        # Add value labels on bars for HIGH severity (most important)
        for i, (h, m) in enumerate(zip(high, moderate)):
            if h > 0:  # Only show if there are HIGH severity interactions
                ax.text(i, h/2, str(h), ha='center', va='center', 
                       fontweight='bold', fontsize=10, color='white')
        
        ax.set_xlabel('Pharmacogenes', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Drug Interactions', fontsize=12, fontweight='bold')
        ax.set_title('Drug-Drug Interaction Severity by Gene\n(Red bars = HIGH risk combinations to avoid)', 
                     fontsize=13, fontweight='bold', pad=15)
        ax.set_xticks(x)
        ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=10)
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=250, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Severity distribution chart saved to: {output_path}")


class GeneDDIVisualizer:
    """Generate gene-specific DDI visualizations for PDF embedding"""
    
    @staticmethod
    def generate_gene_ddi_graph(gene: str, ddis: List[Dict[str, Any]], 
                               output_path: str, max_drugs: int = 15) -> None:
        """Generate a bar chart showing top DDI-involved drugs for a specific gene
        
        Args:
            gene: Gene name
            ddis: List of DDI dictionaries for the gene
            output_path: Path to save the graph PNG file
            max_drugs: Maximum number of top drugs to display
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError as e:
            print(f"⚠️  Warning: matplotlib not available for {gene} graph: {e}")
            return
        
        if not ddis:
            return
        
        # Count interactions per drug
        drug_interaction_count = {}
        drug_severity = {}  # Track highest severity for each drug
        
        severity_weight = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'UNKNOWN': 0.5}
        
        for ddi in ddis:
            drug1 = ddi.get('drug_1', '').strip()
            drug2 = ddi.get('drug_2', '').strip()
            severity = ddi.get('severity', 'UNKNOWN')
            severity_val = severity_weight.get(severity, 0.5)
            
            for drug in [drug1, drug2]:
                if drug:
                    drug_interaction_count[drug] = drug_interaction_count.get(drug, 0) + 1
                    # Track highest severity
                    if drug not in drug_severity:
                        drug_severity[drug] = severity
                    else:
                        if severity_weight.get(severity, 0.5) > severity_weight.get(drug_severity[drug], 0.5):
                            drug_severity[drug] = severity
        
        if not drug_interaction_count:
            return
        
        # Select top drugs
        top_drugs = sorted(drug_interaction_count.keys(), 
                          key=lambda d: drug_interaction_count[d], 
                          reverse=True)[:max_drugs]
        
        # Get counts and colors for top drugs
        counts = [drug_interaction_count[d] for d in top_drugs]
        colors_list = [GeneDDIVisualizer._get_severity_color(drug_severity.get(d, 'UNKNOWN')) 
                      for d in top_drugs]
        
        # Shorten drug names if too long
        display_names = [d[:25] + '...' if len(d) > 25 else d for d in top_drugs]
        
        # Create bar chart
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            bars = ax.barh(range(len(top_drugs)), counts, color=colors_list, edgecolor='black', linewidth=0.5)
            
            ax.set_yticks(range(len(top_drugs)))
            ax.set_yticklabels(display_names, fontsize=10)
            ax.set_xlabel('Number of Interactions', fontsize=11, fontweight='bold')
            ax.set_title(f'{gene} - Top Drugs Involved in DDIs\n({len(ddis)} total interactions)', 
                        fontsize=12, fontweight='bold', pad=15)
            ax.invert_yaxis()
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            
            # Add count labels on bars
            for idx, (bar, count) in enumerate(zip(bars, counts)):
                ax.text(count + max(counts)*0.01, idx, str(count), 
                       va='center', fontsize=9, fontweight='bold')
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='#d62728', edgecolor='black', label='High Severity'),
                Patch(facecolor='#ff7f0e', edgecolor='black', label='Moderate Severity'),
                Patch(facecolor='#ffdb58', edgecolor='black', label='Low Severity'),
                Patch(facecolor='#9467bd', edgecolor='black', label='Unknown Severity')
            ]
            ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Error generating DDI graph for {gene}: {e}")
            return
    
    @staticmethod
    def _get_severity_color(severity: str) -> str:
        """Map severity to color"""
        colors = {
            'HIGH': '#d62728',      # Red
            'MODERATE': '#ff7f0e',  # Orange
            'LOW': '#ffdb58',       # Yellow
            'UNKNOWN': '#9467bd'    # Purple
        }
        return colors.get(severity, '#9467bd')
    
    @staticmethod
    def generate_all_gene_graphs(results: Dict[str, Any], output_prefix: str) -> Dict[str, str]:
        """Generate DDI graphs for all genes
        
        Args:
            results: Analysis results with clinical_decision_support
            output_prefix: Prefix for output files (e.g., 'sample_id')
            
        Returns:
            Dict mapping gene names to their graph file paths
        """
        cds = results.get('clinical_decision_support', {})
        gene_annotations = cds.get('gene_drug_annotations', {})
        
        gene_graph_paths = {}
        
        for gene, gene_data in gene_annotations.items():
            ddis = gene_data.get('drugbank_ddis', [])
            
            if ddis:
                output_path = f"{output_prefix}_ddi_{gene}.png"
                GeneDDIVisualizer.generate_gene_ddi_graph(gene, ddis, output_path, max_drugs=12)
                gene_graph_paths[gene] = output_path
        
        if gene_graph_paths:
            print(f"✓ Generated {len(gene_graph_paths)} gene-specific DDI graphs")
        
        return gene_graph_paths


class GeneDDIPairHeatmap:
    """Generate heatmaps showing which specific drugs interact with which"""
    
    @staticmethod
    def generate_ddi_pair_heatmap(gene: str, ddis: List[Dict[str, Any]], 
                                  output_path: str, top_drugs: int = 20) -> None:
        """Generate a heatmap showing drug-drug interaction pairs
        
        Args:
            gene: Gene name
            ddis: List of DDI dictionaries
            output_path: Path to save PNG
            top_drugs: Maximum drugs to display
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
            import pandas as pd
        except ImportError as e:
            print(f"⚠️  Warning: Required libraries not available for {gene} heatmap: {e}")
            return
        
        if not ddis:
            return
        
        # Get top drugs by interaction frequency
        drug_counts = {}
        for ddi in ddis:
            drug1 = ddi.get('drug_1', '').strip()
            drug2 = ddi.get('drug_2', '').strip()
            if drug1:
                drug_counts[drug1] = drug_counts.get(drug1, 0) + 1
            if drug2:
                drug_counts[drug2] = drug_counts.get(drug2, 0) + 1
        
        if not drug_counts:
            return
        
        # Select top drugs
        top_drug_list = sorted(drug_counts.keys(), 
                              key=lambda d: drug_counts[d], 
                              reverse=True)[:top_drugs]
        
        # Build interaction matrix with severity weighting
        interaction_matrix = {}
        severity_weight = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'UNKNOWN': 0.5}
        
        for drug in top_drug_list:
            interaction_matrix[drug] = {}
            for other_drug in top_drug_list:
                interaction_matrix[drug][other_drug] = 0
        
        # Count interactions between each pair
        for ddi in ddis:
            drug1 = ddi.get('drug_1', '').strip()
            drug2 = ddi.get('drug_2', '').strip()
            severity = ddi.get('severity', 'UNKNOWN')
            weight = severity_weight.get(severity, 0.5)
            
            if drug1 in interaction_matrix and drug2 in interaction_matrix:
                # Since DDIs are bidirectional, add to both directions
                interaction_matrix[drug1][drug2] += weight
                interaction_matrix[drug2][drug1] += weight
        
        # Convert to DataFrame
        df = pd.DataFrame(interaction_matrix).fillna(0)
        
        # Shorten drug names for display
        display_names = []
        for drug in df.index:
            if len(drug) > 20:
                display_names.append(drug[:17] + '...')
            else:
                display_names.append(drug)
        
        # Create heatmap
        try:
            fig, ax = plt.subplots(figsize=(12, 10))
            
            sns.heatmap(df, cmap='RdYlGn_r', annot=True, fmt='.0f',
                       cbar_kws={'label': 'Interaction Severity Score'},
                       linewidths=0.5, linecolor='white', ax=ax,
                       xticklabels=display_names, yticklabels=display_names,
                       vmin=0, vmax=df.values.max())
            
            plt.title(f'{gene} - Drug-Drug Interaction Heatmap\n(Which drug interacts with which)', 
                     fontsize=13, fontweight='bold', pad=15)
            plt.xlabel('Drugs', fontsize=11, fontweight='bold')
            plt.ylabel('Drugs', fontsize=11, fontweight='bold')
            
            # Rotate labels for readability
            plt.xticks(rotation=45, ha='right', fontsize=9)
            plt.yticks(rotation=0, fontsize=9)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error generating heatmap for {gene}: {e}")
            return
    
    @staticmethod
    def generate_all_pair_heatmaps(results: Dict[str, Any], output_prefix: str) -> Dict[str, str]:
        """Generate pair heatmaps for all genes
        
        Args:
            results: Analysis results with clinical_decision_support
            output_prefix: Prefix for output files
            
        Returns:
            Dict mapping gene names to heatmap file paths
        """
        cds = results.get('clinical_decision_support', {})
        gene_annotations = cds.get('gene_drug_annotations', {})
        
        gene_heatmap_paths = {}
        
        for gene, gene_data in gene_annotations.items():
            ddis = gene_data.get('drugbank_ddis', [])
            
            if ddis:
                output_path = f"{output_prefix}_ddi_heatmap_{gene}.png"
                GeneDDIPairHeatmap.generate_ddi_pair_heatmap(gene, ddis, output_path, top_drugs=20)
                gene_heatmap_paths[gene] = output_path
        
        if gene_heatmap_paths:
            print(f"✓ Generated {len(gene_heatmap_paths)} gene-specific DDI pair heatmaps")
        
        return gene_heatmap_paths


if __name__ == "__main__":
    # Example usage
    print("Drug Interaction Integration Module loaded successfully")
    print("Available classes:")
    print("  - DrugInteractionFormatter: Format narratives in NutriGeneEngine style")
    print("  - DrugInteractionExtractor: Extract interactions from analysis results")
    print("  - DrugInteractionReportWriter: Write reports in various formats")
