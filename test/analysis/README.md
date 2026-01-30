# Paper Analysis Code

This folder contains standalone analysis code for the research paper.
These scripts are separate from the main MetaPilot software and can be packaged independently for publication.

## Files

- `PeptideCategoryAnalyzer.java` - Generates peptide distribution statistics for the 4 marker categories
- `run_peptide_analysis.bat` - Batch file to compile and run the analysis

## Usage

### Generate Peptide Distribution Statistics

Run the batch file or compile manually:

```bash
# Compile
javac -d . PeptideCategoryAnalyzer.java

# Run
java analysis.PeptideCategoryAnalyzer
```

The program generates two output files:

1. **`humangut_4category_peptides.csv`** - All peptide distribution statistics (same format as protein distribution)
2. **`humangut_4category_distinct_peptides.csv`** - Distinct peptide counts (peptides found in only ONE genome)

### Memory Management

The distinct peptide analysis is memory-intensive. The program will:
1. First attempt to process all data in a single pass
2. If `OutOfMemoryError` occurs, automatically switch to batch processing by peptide prefix (A-Z)
3. Batch processing reads files multiple times but uses less memory

If you have limited RAM, you can modify the code to use batch processing directly.

## Output Format

CSV with columns:
- Filename
- Total_Peptides
- Anchor, HighResolution, Abundance, Stability (counts)
- Combined (sum of all 4)
- Anchor%, HighRes%, Abundance%, Stability%, Combined% (percentages)

