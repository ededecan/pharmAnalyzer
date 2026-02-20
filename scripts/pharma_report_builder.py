"""
Report generation utilities for pharmacogenomics analyzer.
Centralizes PDF/JSON formatting logic.
"""

import json
from pathlib import Path
from typing import Dict, Any, List
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT


class ReportFormatting:
    """Centralized report formatting constants and utilities"""
    
    # Style definitions
    HEADER_COLOR = colors.HexColor('#1F4E78')
    ALTERNATE_ROW_COLOR = colors.HexColor('#EBF0F7')
    ACCENT_COLOR = colors.HexColor('#2E75B6')
    BORDER_COLOR = colors.HexColor('#C7D3E3')
    TEXT_MUTED = colors.HexColor('#5A6B7A')
    
    # Severity colors (Red for High, Yellow/Orange for Medium, Green for Normal)
    HIGH_PRIORITY_COLOR = colors.HexColor('#C0392B')  # Red
    MEDIUM_PRIORITY_COLOR = colors.HexColor('#F39C12')  # Yellow/Orange
    NORMAL_PRIORITY_COLOR = colors.HexColor('#27AE60')  # Green
    
    # Column widths for different table types
    SUMMARY_TABLE_COLS = [0.8*inch, 1.1*inch, 2.5*inch, 1.5*inch]
    ALLELE_TABLE_COLS = [1.0*inch, 2.0*inch, 2.5*inch]
    
    @staticmethod
    def get_table_style() -> list:
        """Get standard table styling"""
        return [
            ('BACKGROUND', (0, 0), (-1, 0), ReportFormatting.HEADER_COLOR),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('GRID', (0, 0), (-1, -1), 0.5, ReportFormatting.BORDER_COLOR),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, ReportFormatting.ALTERNATE_ROW_COLOR]),
            ('TOPPADDING', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 10),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
            ('RIGHTPADDING', (0, 0), (-1, -1), 8)
        ]


class JSONReportWriter:
    """Handles JSON report generation and output"""
    
    @staticmethod
    def write(analysis_results: Dict[str, Any], output_path: Path) -> None:
        """Write analysis results to JSON file"""
        with open(output_path, 'w') as f:
            json.dump(analysis_results, f, indent=2)
        print(f"\n✓ JSON report saved to: {output_path}")


class PDFReportBuilder:
    """Builds PDF reports with structured content"""
    
    def __init__(self, output_path: str):
        """Initialize PDF builder"""
        self.doc = SimpleDocTemplate(
            output_path,
            pagesize=letter,
            rightMargin=0.6*inch,
            leftMargin=0.6*inch,
            topMargin=0.6*inch,
            bottomMargin=0.6*inch
        )
        self.story = []
        self.styles = getSampleStyleSheet()
        self._setup_custom_styles()
    
    def _setup_custom_styles(self) -> None:
        """Define custom paragraph styles"""
        self.styles.add(ParagraphStyle(
            name='GeneHeader',
            parent=self.styles['Heading2'],
            fontSize=12,
            textColor=ReportFormatting.ACCENT_COLOR,
            spaceAfter=10,
            spaceBefore=14,
            keepWithNext=True
        ))
        
        self.styles.add(ParagraphStyle(
            name='AlleleInfo',
            parent=self.styles['Normal'],
            fontSize=9,
            spaceAfter=6,
            spaceBefore=3,
            leftIndent=12
        ))
        
        self.styles.add(ParagraphStyle(
            name='SectionBreak',
            parent=self.styles['Normal'],
            fontSize=10,
            spaceAfter=8,
            spaceBefore=12,
            keepWithNext=True
        ))

        # Subtle helper styles for metadata and badges
        self.styles.add(ParagraphStyle(
            name='Muted',
            parent=self.styles['Normal'],
            fontSize=8.5,
            textColor=ReportFormatting.TEXT_MUTED,
            spaceAfter=4
        ))

        self.styles.add(ParagraphStyle(
            name='Small',
            parent=self.styles['Normal'],
            fontSize=8.5,
            spaceAfter=4
        ))

        self.styles.add(ParagraphStyle(
            name='Badge',
            parent=self.styles['Normal'],
            fontSize=8,
            textColor=colors.whitesmoke,
            backColor=ReportFormatting.ACCENT_COLOR,
            leftIndent=0,
            rightIndent=0,
            leading=10,
            spaceAfter=4,
            spaceBefore=0
        ))

        self.styles.add(ParagraphStyle(
            name='EvidenceNote',
            parent=self.styles['Normal'],
            fontSize=8.7,
            leftIndent=18,
            textColor=ReportFormatting.TEXT_MUTED,
            spaceAfter=6,
            spaceBefore=0,
            leading=11
        ))
        
        self.styles.add(ParagraphStyle(
            name='VariantNotation',
            parent=self.styles['Normal'],
            fontSize=9,
            fontName='Courier',
            leftIndent=12,
            spaceAfter=4
        ))
        
        # Therapeutic area and drug interaction styles
        self.styles.add(ParagraphStyle(
            name='TherapeuticArea',
            parent=self.styles['Normal'],
            fontSize=10,
            textColor=colors.HexColor('#1F4E78'),
            spaceAfter=6,
            spaceBefore=4,
            leftIndent=6,
            fontName='Helvetica-Bold'
        ))
        
        self.styles.add(ParagraphStyle(
            name='DrugName',
            parent=self.styles['Normal'],
            fontSize=9,
            textColor=ReportFormatting.ACCENT_COLOR,
            spaceAfter=4,
            fontName='Helvetica-Bold'
        ))
        
        self.styles.add(ParagraphStyle(
            name='DiplotypeEvidence',
            parent=self.styles['Normal'],
            fontSize=8.8,
            spaceAfter=8,
            spaceBefore=4,
            leftIndent=12,
            leading=11
        ))
        
        self.styles.add(ParagraphStyle(
            name='ClinicalEvidenceSeq',
            parent=self.styles['Normal'],
            fontSize=8.5,
            spaceAfter=6,
            leftIndent=12,
            leading=10
        ))
    
    def add_title(self, title: str) -> None:
        """Add report title"""
        self.story.append(Paragraph(f"<b>{title}</b>", self.styles['Heading1']))
        self.story.append(Spacer(1, 0.2*inch))
    
    def add_section_header(self, header: str) -> None:
        """Add section header"""
        self.story.append(Paragraph(f"<b>{header}</b>", self.styles['Heading2']))
        self.story.append(Spacer(1, 0.1*inch))
    
    def add_subsection_header(self, header: str) -> None:
        """Add subsection header"""
        self.story.append(Paragraph(f"<b>{header}</b>", self.styles['Heading3']))
        self.story.append(Spacer(1, 0.08*inch))
    
    def add_table(self, data: list, col_widths: list) -> None:
        """Add formatted table"""
        table = Table(data, colWidths=col_widths)
        table.setStyle(TableStyle(ReportFormatting.get_table_style()))
        self.story.append(table)
    
    def add_paragraph(self, text: str, style_name: str = 'Normal') -> None:
        """Add paragraph text"""
        self.story.append(Paragraph(text, self.styles[style_name]))
    
    def add_image(self, image_path: str, width: float = 7.0, height: float = 5.5) -> None:
        """Add image to PDF"""
        try:
            from reportlab.platypus import Image as RLImage
            img = RLImage(image_path, width=width*inch, height=height*inch)
            self.story.append(img)
        except Exception as e:
            print(f"Warning: Could not embed image {image_path}: {e}")
    
    def add_spacer(self, height: float = 0.1) -> None:
        """Add vertical spacer"""
        self.story.append(Spacer(1, height*inch))
    
    def get_story_height(self) -> int:
        """Get approximate height of current story content"""
        return len(self.story) * 12  # Rough estimate in points
    
    def add_page_break(self) -> None:
        """Add page break"""
        self.story.append(PageBreak())
    
    def build(self) -> None:
        """Build the PDF document"""
        self.doc.build(self.story)
    
    def get_story(self) -> list:
        """Get story content for manual manipulation"""
        return self.story


class ConsolidatedReportHelper:
    """Helper class for building consolidated reports with therapeutic areas and heatmaps"""
    
    @staticmethod
    def get_therapeutic_area_badge(drug_name: str) -> str:
        """Get therapeutic area for a drug with color coding"""
        try:
            from atc_classification import ATCClassifier
            area = ATCClassifier.get_therapeutic_area(drug_name)
            return f"<font color='#2E75B6'>[{area.value}]</font>"
        except:
            return ""
    
    @staticmethod
    def format_drug_with_area(drug_name: str, therapeutic_area: str = None) -> str:
        """Format drug name with therapeutic area"""
        badge = ConsolidatedReportHelper.get_therapeutic_area_badge(drug_name)
        if badge:
            return f"{badge} <b>{drug_name}</b>"
        return f"<b>{drug_name}</b>"
    
    @staticmethod
    def get_drugs_by_therapeutic_area(drug_list: List[str]) -> Dict[str, List[str]]:
        """Organize drugs by therapeutic area"""
        try:
            from atc_classification import ATCClassifier
            return ATCClassifier.organize_drugs_by_area(drug_list)
        except:
            return {drug: [drug] for drug in drug_list}


class SummaryTableBuilder:
    """Builds summary tables from analysis results"""
    
    @staticmethod
    def build_summary_table(analysis_results: Dict[str, Any]) -> tuple:
        """Build summary table data and return with column widths"""
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import Paragraph
        
        styles = getSampleStyleSheet()
        header_style = styles['Normal'].clone('HeaderStyle')
        header_style.fontSize = 10
        header_style.textColor = colors.whitesmoke
        header_style.alignment = 0  # LEFT
        
        cell_style = styles['Normal'].clone('CellStyle')
        cell_style.fontSize = 9.5
        cell_style.leading = 12
        
        # Header row with Paragraph objects
        data = [[
            Paragraph('<b>Gene</b>', header_style),
            Paragraph('<b>Diplotype</b>', header_style),
            Paragraph('<b>Phenotype</b>', header_style),
            Paragraph('<b>EHR Priority</b>', header_style)
        ]]
        
        for gene in sorted(analysis_results['variants_by_gene'].keys()):
            gene_data = analysis_results['variants_by_gene'][gene]
            data.append([
                Paragraph(f"<b>{gene}</b>", cell_style),
                Paragraph(str(gene_data.get('diplotype') or 'N/A'), cell_style),
                Paragraph(str(gene_data.get('coded_summary') or 'N/A'), cell_style),
                Paragraph(str(gene_data.get('ehr_priority') or 'none'), cell_style)
            ])
        
        return data, ReportFormatting.SUMMARY_TABLE_COLS
