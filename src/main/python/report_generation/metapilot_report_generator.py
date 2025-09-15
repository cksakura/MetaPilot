from .base_generator import BaseReportGenerator
from .id_summary_generator import IDSummaryGenerator
from .peptides_summary_generator import PeptidesSummaryGenerator
from .proteins_summary_generator import ProteinsSummaryGenerator
from .taxonomy_summary_generator import TaxonomySummaryGenerator
from .function_summary_generator import FunctionSummaryGenerator

class MetaPilotReportGenerator(BaseReportGenerator):
    """Enhanced report generator for MetaPilot TSV file formats with metadata support."""
    
    def __init__(self, output_dir="reports"):
        super().__init__(output_dir)
        
        # Initialize specialized generators
        self.id_generator = IDSummaryGenerator(output_dir)
        self.peptides_generator = PeptidesSummaryGenerator(output_dir)
        self.proteins_generator = ProteinsSummaryGenerator(output_dir)
        self.taxonomy_generator = TaxonomySummaryGenerator(output_dir)
        self.function_generator = FunctionSummaryGenerator(output_dir)
    
    def create_id_summary(self, summary_df, metadata=None):
        """Generate ID summary report with metadata support."""
        return self.id_generator.create_id_summary(summary_df, metadata)
    
    def create_peptides_summary(self, peptides_df, metadata=None):
        """Enhanced peptides summary with metadata support."""
        return self.peptides_generator.create_peptides_summary(peptides_df, metadata)
    
    def create_proteins_summary(self, proteins_df, metadata=None):
        """Enhanced proteins summary with advanced dimensionality reduction."""
        return self.proteins_generator.create_proteins_summary(proteins_df, metadata)
    
    def create_taxonomy_summary(self, taxa_df, metadata=None):
        """Enhanced taxonomy summary with rank-specific analysis for all taxonomic ranks."""
        return self.taxonomy_generator.create_taxonomy_summary(taxa_df, metadata)
    
    def create_function_summary(self, func_df, metadata=None):
        """Enhanced function summary with category analysis and dimensionality reduction."""
        return self.function_generator.create_function_summary(func_df, metadata)