> Historical archive. This prompt is retained for project history. For the current workflow strategy, use [../user/workflows.md](../user/workflows.md).

Create a polished scientific workflow diagram for a bioinformatics software paper and presentation. 

Style:
- clean publication-quality vector infographic
- professional scientific figure, not cartoon, not pixel art
- white background, subtle light gray panel shading
- crisp black text, minimal gradients, thin clean arrows
- modern journal figure style, compact but readable
- balanced layout, high visual hierarchy, no clutter
- all labels large and legible
- rectangular rounded boxes, consistent alignment
- use only 3 pathway colors with soft outlines:
  - blue for DDA
  - green for DIA Marker
  - orange for DIA Seed
- shared core steps in neutral gray
- include a small clean legend in the top right

Overall layout:
A wide horizontal workflow with 5 numbered phases from top to bottom or left to right, clearly separated as large horizontal bands:
1. Input
2. Initial evidence search
3. Genome ranking & selection
4. Refined search-space construction
5. Final quantification & export

Phase 1: INPUT
Show:
- Raw MS data (DDA / DIA)
Keep this simple and compact.

Phase 2: INITIAL EVIDENCE SEARCH (small search space)
Split into DDA and DIA branches:
- DDA branch (blue): “High-abundance / Marker FASTA” → “Initial PSMs & peptides”
- DIA Marker branch (green): “Marker spectral library” → “Initial precursor evidence”
- DIA Seed branch (orange): “User-selected libraries” → “Initial precursor evidence”
Then merge these into one shared box:
- “Evidence from marker-associated peptides / proteins”

Phase 3: GENOME RANKING & SELECTION
This should be the visual center of the figure.
Show four connected conceptual modules:
1. “Evidence extraction”
   subtitle: peptides + proteins + MAG metadata
2. “Seed genomes (optional)”
   subtitle: privileged genomes + taxonomic representatives
3. “Greedy genome ordering”
   subtitle: maximize incremental peptide coverage
   include a simple clean curve icon
4. “Adaptive knee-stop selection”
   subtitle: Focus vs Comprehensive
   show one small curve with two stop points:
   - Focus stop = earlier, stricter
   - Comprehensive stop = later, broader
Output box at bottom:
- “Selected genomes”

Phase 4: REFINED SEARCH-SPACE CONSTRUCTION
From “Selected genomes”, split into:
- DDA branch (blue): “Sample-specific FASTA (selected genomes + optional host)”
- DIA branches (green/orange): “Refined microbiome library”
Show that both DIA Marker and DIA Seed use the refined microbiome library after genome selection.
Also visually indicate that the Seed workflow begins with user libraries but converges to the same refined library logic.

Phase 5: FINAL QUANTIFICATION & EXPORT
Show:
- “Refined quantifying search”
- arrow to “Final IDs & intensities”
- arrow to “Unified taxonomy and function export”
Final export box subtitle:
- “Genome, species and functions”

Include one small side note or callout for Deep workflow:
- “Deep mode: iterative rounds of temporary refined library search before final quantification”
Place this as a neat inset near the DIA branch, not dominant.
Use a small loop arrow icon to indicate repeated rounds.
Do not let Deep dominate the full diagram.

Important content rules:
- emphasize that all workflows share one common backbone
- show that only the initial evidence search space differs
- DDA uses high-abundance / marker FASTA
- DIA Marker uses Marker spectral library
- DIA Seed uses user-selected libraries
- genome ranking and selection is shared logic
- refined search space is built from selected genomes
- Focus and Comprehensive differ only by adaptive stopping stringency
- Deep adds iterative refinement rounds before the final quantifying search

Design constraints:
- avoid overcrowding
- avoid too many file names
- avoid long paragraphs
- avoid cartoon scientists, avatars, trophies, or decorative mascots
- avoid pixel art
- avoid unnecessary icons
- keep the figure suitable for both a methods paper and a conference talk
- make it look like a Nature Methods / Bioinformatics style summary figure

Text to include exactly where helpful:
- “Initial evidence search (small search space)”
- “Genome ranking & selection”
- “Refined search-space construction”
- “Refined quantifying search”
- “Unified taxonomy export”
- “Focus: earlier stop”
- “Comprehensive: later stop”
- “Deep: iterative rounds before final search”

Output:
- single high-resolution landscape scientific workflow figure
- publication-ready, presentation-ready, visually clean, readable at slide scale
