import pandas as pd
from report_generator import ReportGenerator

def main():
    # Create sample data
    data = pd.DataFrame({
        'Protein': ['P1', 'P2', 'P3', 'P4', 'P5'],
        'Intensity': [1000, 2000, 1500, 3000, 2500],
        'Peptides': [5, 8, 6, 10, 7],
        'Coverage': [0.3, 0.5, 0.4, 0.6, 0.45]
    })
    
    # Initialize report generator
    generator = ReportGenerator(output_dir="example_reports")
    
    # Define plot configurations
    plots = [
        {
            "type": "correlation",
            "output": "protein_correlation.html"
        },
        {
            "type": "distribution",
            "column": "Intensity",
            "output": "intensity_distribution.html"
        },
        {
            "type": "distribution",
            "column": "Coverage",
            "output": "coverage_distribution.html"
        }
    ]
    
    # Generate report
    report = generator.generate_report(
        data,
        "Protein Analysis Report",
        plots
    )
    
    print(f"Report generated: {report['name']}")
    print("\nGenerated plots:")
    for plot in report['plots']:
        print(f"- {plot['type']} plot: {plot['path']}")

if __name__ == "__main__":
    main() 