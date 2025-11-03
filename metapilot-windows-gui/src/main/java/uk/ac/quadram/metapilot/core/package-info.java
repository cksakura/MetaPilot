/**
 * MetaPilot Core Package - Fundamental Bioinformatics Components
 * 
 * <h2>Overview</h2>
 * The core package provides fundamental bioinformatics components and utilities 
 * that form the foundation of the MetaPilot metaproteomics analysis platform.
 * These components handle the basic building blocks of protein and peptide analysis.
 * 
 * <h2>Package Structure</h2>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.aminoacid} - Amino Acid Management</h3>
 * <ul>
 *   <li>{@link uk.ac.quadram.metapilot.core.aminoacid.Aminoacids} - Comprehensive amino acid properties and modifications</li>
 *   <li>{@link uk.ac.quadram.metapilot.core.aminoacid.Aminoacid} - Individual amino acid representation</li>
 *   <li>{@link uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty} - Physical and chemical constants</li>
 *   <li>{@link uk.ac.quadram.metapilot.core.aminoacid.AminoacidFragment} - Peptide fragmentation calculations</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.enzyme} - Proteolytic Enzyme Support</h3>
 * <ul>
 *   <li>{@link uk.ac.quadram.metapilot.core.enzyme.Enzyme} - Protease specificity and digestion simulation</li>
 *   <li>Support for trypsin, LysC, ArgC, chymotrypsin, and other common proteases</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.function} - Functional Annotation</h3>
 * <ul>
 *   <li>Protein function prediction and classification</li>
 *   <li>COG, KEGG, and enzyme commission number support</li>
 *   <li>Functional database integration</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.math} - Mathematical Utilities</h3>
 * <ul>
 *   <li>Statistical analysis functions</li>
 *   <li>Numerical algorithms for proteomics</li>
 *   <li>Distribution calculations</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.MGYG} - MGnify Integration</h3>
 * <ul>
 *   <li>MGnify database access and management</li>
 *   <li>Metagenomic data integration</li>
 *   <li>FTP download utilities</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.ml} - Machine Learning</h3>
 * <ul>
 *   <li>Machine learning models for proteomics</li>
 *   <li>Feature extraction and selection</li>
 *   <li>Prediction algorithms</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.mod} - Modification Management</h3>
 * <ul>
 *   <li>{@link uk.ac.quadram.metapilot.core.mod.UnimodReader} - Unimod database integration</li>
 *   <li>Post-translational modification support</li>
 *   <li>Chemical modification definitions</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.model} - Data Models</h3>
 * <ul>
 *   <li>Core data structures for proteomics objects</li>
 *   <li>Protein, peptide, and spectrum representations</li>
 *   <li>Result container classes</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.prodb} - Protein Database</h3>
 * <ul>
 *   <li>Protein database management</li>
 *   <li>FASTA file processing</li>
 *   <li>Sequence indexing and search</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.prosit} - Prosit Integration</h3>
 * <ul>
 *   <li>Deep learning-based spectrum prediction</li>
 *   <li>Prosit model integration</li>
 *   <li>Spectral library generation</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.taxonomy} - Taxonomic Classification</h3>
 * <ul>
 *   <li>NCBI taxonomy integration</li>
 *   <li>Taxonomic tree navigation</li>
 *   <li>Species and strain identification</li>
 * </ul>
 * 
 * <h3>{@link uk.ac.quadram.metapilot.core.tools} - Utility Tools</h3>
 * <ul>
 *   <li>{@link uk.ac.quadram.metapilot.core.tools.MwCalculator} - Molecular weight calculations</li>
 *   <li>File processing utilities</li>
 *   <li>Conversion tools</li>
 * </ul>
 * 
 * <h2>Design Principles</h2>
 * <ul>
 *   <li><strong>Modularity:</strong> Self-contained components with clear interfaces</li>
 *   <li><strong>Reusability:</strong> Generic implementations usable across different workflows</li>
 *   <li><strong>Performance:</strong> Optimized algorithms for large-scale data processing</li>
 *   <li><strong>Extensibility:</strong> Easy to extend with new functionality</li>
 * </ul>
 * 
 * <h2>Key Features</h2>
 * <ul>
 *   <li><strong>Mass Accuracy:</strong> High-precision mass calculations for both monoisotopic and average masses</li>
 *   <li><strong>Modification Support:</strong> Comprehensive support for chemical modifications</li>
 *   <li><strong>Multi-organism:</strong> Support for diverse microbial and metagenomic datasets</li>
 *   <li><strong>Standards Compliance:</strong> Adherence to proteomics standards and ontologies</li>
 * </ul>
 * 
 * <h2>Dependencies</h2>
 * <p>The core package is designed to minimize external dependencies while providing
 * rich functionality. Key dependencies include:</p>
 * <ul>
 *   <li>Apache Commons for utilities</li>
 *   <li>DOM4J for XML processing</li>
 *   <li>Standard Java libraries</li>
 * </ul>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 */
package uk.ac.quadram.metapilot.core;
