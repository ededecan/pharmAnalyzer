#!/usr/bin/env python3
"""
Gene Prioritization System for Pharmacogenomics Analysis
Based on clinical importance and PharmGKB VIP status
"""

import json
from typing import Dict, List, Tuple

class GenePrioritizer:
    """Classifies pharmacogenes by clinical priority levels"""
    
    def __init__(self):
        self.priority_classes = {
            "CRITICAL_SAFETY": {
                "level": 1,
                "description": "Life-threatening toxicity prevention",
                "genes": ["DPYD", "TPMT", "NUDT15", "G6PD", "UGT1A1"],
                "clinical_action": "MANDATORY_TESTING",
                "color_code": "RED"
            },
            "MAJOR_METABOLIZERS": {
                "level": 2, 
                "description": "Drug metabolism >40% of clinical drugs",
                "genes": ["CYP2C9", "CYP2C19", "CYP2D6"],
                "clinical_action": "STRONGLY_RECOMMENDED", 
                "color_code": "ORANGE"
            },
            "DRUG_TRANSPORTERS": {
                "level": 3,
                "description": "Drug uptake/efflux controllers",
                "genes": ["SLCO1B1", "ABCB1", "ABCG2"],
                "clinical_action": "RECOMMENDED",
                "color_code": "YELLOW"
            },
            "SPECIALIZED_THERAPY": {
                "level": 4,
                "description": "Specific drug classes or conditions", 
                "genes": ["VKORC1", "CYP2B6", "IFNL3", "CYP3A5", "CYP4F2"],
                "clinical_action": "CONDITIONAL",
                "color_code": "GREEN"
            },
            "CLOTTING_FACTORS": {
                "level": 5,
                "description": "Thrombophilia & anticoagulation guidance",
                "genes": ["F2", "F5"],
                "clinical_action": "RISK_ASSESSMENT",
                "color_code": "BLUE"
            },
            "EXTENDED_PANEL": {
                "level": 6,
                "description": "Additional CYP enzymes for comprehensive analysis",
                "genes": ["CYP1A1", "CYP1A2", "CYP2C8", "CYP3A4", "ALDH2", "BCHE", "COMT"],
                "clinical_action": "RESEARCH_INFORMATIVE",
                "color_code": "GRAY"
            }
        }
    
    def get_gene_priority(self, gene: str) -> Dict:
        """Get priority classification for a specific gene"""
        for category, info in self.priority_classes.items():
            if gene in info["genes"]:
                return {
                    "gene": gene,
                    "category": category,
                    "level": info["level"],
                    "description": info["description"],
                    "clinical_action": info["clinical_action"],
                    "color_code": info["color_code"]
                }
        return {"gene": gene, "category": "UNKNOWN", "level": 99}
    
    def prioritize_gene_list(self, genes: List[str]) -> List[Dict]:
        """Sort genes by clinical priority"""
        prioritized = []
        for gene in genes:
            priority_info = self.get_gene_priority(gene)
            prioritized.append(priority_info)
        
        # Sort by priority level
        return sorted(prioritized, key=lambda x: x.get("level", 99))
    
    def get_critical_genes(self) -> List[str]:
        """Get genes requiring mandatory testing"""
        return self.priority_classes["CRITICAL_SAFETY"]["genes"]
    
    def generate_priority_report(self, genes: List[str]) -> Dict:
        """Generate comprehensive priority analysis"""
        prioritized = self.prioritize_gene_list(genes)
        
        report = {
            "total_genes": len(genes),
            "priority_breakdown": {},
            "clinical_recommendations": [],
            "gene_details": prioritized
        }
        
        # Count genes by category
        for category in self.priority_classes.keys():
            category_genes = [g for g in prioritized if g.get("category") == category]
            report["priority_breakdown"][category] = {
                "count": len(category_genes),
                "genes": [g["gene"] for g in category_genes]
            }
        
        # Generate recommendations 
        critical_count = len(report["priority_breakdown"].get("CRITICAL_SAFETY", {}).get("genes", []))
        major_count = len(report["priority_breakdown"].get("MAJOR_METABOLIZERS", {}).get("genes", []))
        
        if critical_count > 0:
            report["clinical_recommendations"].append(
                f"PRIORITY 1: Test {critical_count} critical safety genes to prevent life-threatening toxicity"
            )
        if major_count > 0:
            report["clinical_recommendations"].append(
                f"PRIORITY 2: Test {major_count} major metabolizers covering >40% of clinical drugs"
            )
            
        return report

if __name__ == "__main__":
    # Example usage
    prioritizer = GenePrioritizer()
    
    # Test with our gene list
    genes = [
        "ABCB1", "ABCG2", "ALDH2", "BCHE", "COMT", "CYP1A1", "CYP1A2", 
        "CYP2B6", "CYP2C19", "CYP2C8", "CYP2C9", "CYP2D6", "CYP3A4", 
        "CYP3A5", "CYP4F2", "DPYD", "F2", "F5", "G6PD", "IFNL3", 
        "NUDT15", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"
    ]
    
    report = prioritizer.generate_priority_report(genes)
    print(json.dumps(report, indent=2))