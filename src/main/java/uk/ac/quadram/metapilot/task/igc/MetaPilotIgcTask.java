package uk.ac.quadram.metapilot.task.igc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.math.BigInteger;
import java.util.concurrent.ExecutionException;

import uk.ac.quadram.metapilot.core.taxonomy.Taxon;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabase;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyDatabaseSQL;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPeptide;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.MetaPilotTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportCopyTask;
import uk.ac.quadram.metapilot.task.MetaPythonReportTask;
import uk.ac.quadram.metapilot.task.io.MetaAlgorithm;
import uk.ac.quadram.metapilot.task.io.MetaBiomJsonHandler;
import uk.ac.quadram.metapilot.task.io.MetaTreeHandler;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptideXMLWriter;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinAnnoIgc;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.mdb.unipept.UnipeptRequester;
import uk.ac.quadram.metapilot.mdb.unipept.UnipProResult;
import uk.ac.quadram.metapilot.mdb.unipept.UnipProInfo;
import uk.ac.quadram.metapilot.mdb.unipept.UnipLCAResult;

/**
 * Abstract base class for DDA IGC (Integrated Gene Catalog) workflow tasks.
 * 
 * This class provides common functionality and fields shared across all DDA IGC workflows,
 * including pFind, MaxQuant, Sage, and FragPipe implementations. It handles common operations
 * such as result directory creation, report generation, and taxonomic analysis.
 */
public abstract class MetaPilotIgcTask extends MetaPilotTask {

    // Common fields for all IGC workflows
    protected String[] rawFiles;
    protected String[] rawFileNames;
    protected String[] quanExpNames;
    protected File spectraDirFile;
    protected File[] separateResultFolders;
    protected int threadCount;
    protected double psmCoverThreshold = 0.95;

    protected File rawDir;
    protected File quan_pep_file;
    protected File quan_pro_file;
    protected File final_pep_txt;
    protected File final_pro_txt;
    protected HashMap<String, String> expNameMap;
    protected File summaryMetaFile;
    protected File combinedFile;
    
    protected MetaPeptide[] metaPeptides;

    public MetaPilotIgcTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
            ProgressBar secondBar) {
        super(metaPar, mrp, progressBar, secondBar);

        // Report templates will be copied in background during call() method
        // This prevents UI freeze during task construction
        LOGGER.info(getTaskName() + ": Constructor completed, report templates will be copied during execution");
    }
    
    /**
     * Export report including peptide/protein reports, taxonomy analysis, and functional annotation
     * This method follows the same pattern as DDAIgcPFindTask but uses Sage readers
     */
    protected boolean exportReport() {
        updateMessage(getTaskName() + ": exporting report");

        LOGGER.info(getTaskName() + ": exporting report started");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting report started");

        // Update progress before report generation
        updateProgress(75, 100);
        updateMessage("Starting report generation...");

        File reportDir = new File(this.metaPar.getResult(), "report");
        if (!reportDir.exists()) {
            reportDir.mkdirs();
        }

        // Copy report template
        MetaPythonReportCopyTask reportCopyTask = new MetaPythonReportCopyTask(new File(metaPar.getResult()));
        Thread reportCopyThread = new Thread(reportCopyTask);
        reportCopyThread.start();
        try {
            reportCopyThread.join();
            if (!reportCopyTask.get()) {
                LOGGER.error(getTaskName() + ": error in copying report template");
                System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template");
                return false;
            }
        } catch (InterruptedException | ExecutionException e) {
            LOGGER.error(getTaskName() + ": error in copying report template", e);
            System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in copying report template: " + e.getMessage());
            return false;
        }

        metaPar.setReportDir(reportDir);

        File reportHtmlDir = new File(reportDir + "\\reports\\rmd_html");
        if (!reportHtmlDir.exists()) {
            reportHtmlDir.mkdirs();
        }

        // Export metadata
        File summaryMetaFile = new File(this.metaPar.getResult(), "metainfo.tsv");
        try {
            this.metaPar.getMetadata().exportMetadata(summaryMetaFile);
        } catch (FileNotFoundException e) {
            LOGGER.error("Error exporting metadata", e);
        }

        MetaPythonReportTask task = new MetaPythonReportTask(metaPar.getThreadCount());

        // Export summary
        boolean summaryFinished = this.exportSummary(task);
        updateProgress(80, 100);

        // Generate final peptide and protein files
        String[] rawFiles = this.metaPar.getMetadata().getRawFiles();
        String[] quanExpNames = new String[rawFiles.length];
        for (int i = 0; i < rawFiles.length; i++) {
            quanExpNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1, rawFiles[i].lastIndexOf("."));
        }

        // Use the common processPeptideListAndGenerateReports method
        // This method will handle protein inference and report generation
        // TODO: Build file mappings and call processPeptideListAndGenerateReports

        // For now, keep the existing exportReport logic
        // Generate final peptide report
        File final_pep_txt = new File(this.metaPar.getDbSearchResultFile(), "final_peptides.tsv");
        if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
            LOGGER.info(getTaskName() + ": exporting peptide report started");
            System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

            try {
                
                try (PrintWriter writer = new PrintWriter(final_pep_txt)) {
                    // Write header
                    StringBuilder titlesb = new StringBuilder();
                    titlesb.append("Sequence").append("\t");
                    titlesb.append("Base sequence").append("\t");
                    titlesb.append("Length").append("\t");
                    titlesb.append("Missed cleavages").append("\t");
                    titlesb.append("Mass").append("\t");
                    titlesb.append("Proteins").append("\t");
                    titlesb.append("Charges").append("\t");
                    titlesb.append("Score").append("\t");
                    titlesb.append("Final score").append("\t");
                    titlesb.append("Q-value").append("\t");

                    for (int i = 0; i < quanExpNames.length; i++) {
                        titlesb.append("Identification type " + quanExpNames[i]).append("\t");
                    }

                    titlesb.append("Intensity").append("\t");

                    for (int i = 0; i < quanExpNames.length; i++) {
                        titlesb.append("Intensity " + quanExpNames[i]).append("\t");
                    }

                    titlesb.append("Reverse").append("\t");
                    titlesb.append("Potential contaminant").append("\t");
                    titlesb.append("id").append("\t");
                    titlesb.append("Protein group IDs").append("\t");
                    titlesb.append("MS/MS Count");

                    writer.println(titlesb);

                    // Write peptide data
                    for (int i = 0; i < metaPeptides.length; i++) {
                        MetaPeptide pep = metaPeptides[i];
                        StringBuilder sb = new StringBuilder();
                        sb.append(pep.getSequence()).append("\t");
                        sb.append(pep.getSequence()).append("\t");
                        sb.append(pep.getLength()).append("\t");
                        sb.append(pep.getMissCleave()).append("\t");
                        sb.append(pep.getMass()).append("\t");

                        String[] proteins = pep.getProteins();
                        if (proteins.length > 0) {
                            for (String pro : proteins) {
                                sb.append(pro).append(";");
                            }
                            sb.deleteCharAt(sb.length() - 1);
                        }
                        sb.append("\t");

                        int[] charge = pep.getCharges();
                        StringBuilder chargesb = new StringBuilder();
                        for (int j = 0; j < charge.length; j++) {
                            if (charge[j] == 1) {
                                chargesb.append(j + 1).append(";");
                            }
                        }
                        if (chargesb.length() > 1) {
                            sb.append(chargesb.subSequence(0, chargesb.length() - 1)).append("\t");
                        } else {
                            sb.append("2").append("\t");
                        }

                        sb.append(pep.getScore()).append("\t");
                        sb.append(pep.getScore()).append("\t");
                        sb.append(pep.getPEP()).append("\t");

                        // Add identification types and intensities
                        double[] intensity = pep.getIntensity();
                        int[] idenType = pep.getIdenType();
                        
                        for (int j = 0; j < idenType.length; j++) {
                            if (idenType[j] == 0) {
                                sb.append("By MS/MS").append("\t");
                            } else {
                                sb.append("\t");
                            }
                        }

                        double totalIntensity = 0;
                        for (int j = 0; j < intensity.length; j++) {
                            totalIntensity += intensity[j];
                        }
                        sb.append((int) totalIntensity).append("\t");

                        for (int j = 0; j < intensity.length; j++) {
                            sb.append((int) intensity[j]).append("\t");
                        }

                        sb.append(pep.isUse() ? "" : "+").append("\t");
                        sb.append("").append("\t");
                        sb.append(i).append("\t");
                        sb.append("").append("\t");
                        sb.append(pep.getTotalMS2Count());

                        writer.println(sb);
                    }
                }

                LOGGER.info(getTaskName() + ": peptide report has been exported to " + final_pep_txt);
                System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": peptide report has been exported to " + final_pep_txt);

            } catch (Exception e) {
                LOGGER.error(getTaskName() + ": error in reading peptides from Sage results", e);
                System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading peptides from Sage results");
                return false;
            }
        }

        updateProgress(85, 100);

        // Generate final protein report (simplified for Sage)
        File final_pro_txt = new File(this.metaPar.getResult(), "final_proteins.tsv");
        if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
            LOGGER.info(getTaskName() + ": exporting protein report started");
            System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting protein report started");

            try {
                // For Sage, we need to create a simplified protein report from peptide data
                // This is a simplified version - in practice, you'd read the actual Sage protein results
                try (PrintWriter writer = new PrintWriter(final_pro_txt)) {
                    // Write header
                    StringBuilder titlesb = new StringBuilder();
                    titlesb.append("Protein IDs").append("\t");
                    titlesb.append("Majority protein IDs").append("\t");
                    titlesb.append("Peptide counts (all)").append("\t");
                    titlesb.append("Peptide counts (razor)").append("\t");
                    titlesb.append("Number of proteins").append("\t");
                    titlesb.append("Peptides").append("\t");
                    titlesb.append("Score").append("\t");
                    titlesb.append("Intensity").append("\t");

                    for (int i = 0; i < quanExpNames.length; i++) {
                        titlesb.append("Intensity " + quanExpNames[i]).append("\t");
                    }

                    titlesb.append("Reverse").append("\t");
                    titlesb.append("Potential contaminant").append("\t");
                    titlesb.append("id").append("\t");
                    titlesb.append("Peptide IDs").append("\t");
                    titlesb.append("Peptide is razor");

                    writer.println(titlesb);

                    // Create simplified protein data from peptide data
                    // This is a placeholder - in practice, you'd read actual Sage protein results
                    writer.println("Protein_1\tProtein_1\t1\t1\t1\t1\t0.0\t0.0\t" + 
                        String.join("\t", new String[quanExpNames.length]) + "\t\t\t0\t\t");
                }

                LOGGER.info(getTaskName() + ": protein report has been exported to " + final_pro_txt);
                System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": protein report has been exported to " + final_pro_txt);

            } catch (Exception e) {
                LOGGER.error(getTaskName() + ": error in reading proteins from Sage results", e);
                System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in reading proteins from Sage results");
                return false;
            }
        }

        updateProgress(90, 100);

        // Taxonomy analysis (same as DDAIgcPFindTask)
        LOGGER.info(getTaskName() + ": taxonomy analysis started");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": taxonomy analysis started");

        File taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
        if (!taxFile.exists()) {
            taxFile.mkdir();
        }

        // Read peptide intensity data for taxonomy analysis
        HashMap<String, double[]> pepIntensityMap = new HashMap<String, double[]>();
        String[] expNames = null;
        try (BufferedReader quanPepReader = new BufferedReader(new FileReader(final_pep_txt))) {
            String line = quanPepReader.readLine(); // Skip header
            String[] title = line.split("\t");
            int seqId = -1;
            ArrayList<String> list = new ArrayList<String>();
            for (int i = 0; i < title.length; i++) {
                if (title[i].equals("Sequence")) {
                    seqId = i;
                } else if (title[i].startsWith("Intensity ")) {
                    list.add(title[i]);
                }
            }

            expNames = list.toArray(new String[list.size()]);
            int[] intenId = new int[list.size()];
            for (int i = 0; i < expNames.length; i++) {
                for (int j = 0; j < title.length; j++) {
                    if (expNames[i].equals(title[j])) {
                        intenId[i] = j;
                    }
                }
            }

            while ((line = quanPepReader.readLine()) != null) {
                String[] cs = line.split("\t");
                double[] intensity = new double[expNames.length];
                for (int i = 0; i < intensity.length; i++) {
                    intensity[i] = Double.parseDouble(cs[intenId[i]]);
                }
                pepIntensityMap.put(cs[seqId], intensity);
            }
        } catch (IOException e) {
            LOGGER.error("Error in reading peptides from " + final_pep_txt, e);
            return false;
        }

        // Execute Unipept analysis (same as DDAIgcPFindTask)
        boolean executeUnipept = true; // For now, always execute
        MetaPep2TaxaPar mptp = MetaPep2TaxaPar.getDefault();
        File taxResultFile = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".taxonomy.xml");
        
        if (executeUnipept) {
            // Get unique peptide sequences
            HashSet<String> sequenceSet = new HashSet<String>();
            try (BufferedReader reader = new BufferedReader(new FileReader(final_pep_txt))) {
                String line = reader.readLine(); // Skip header
                int seqId = -1;
                String[] title = line.split("\t");
                for (int i = 0; i < title.length; i++) {
                    if (title[i].equals("Sequence")) {
                        seqId = i;
                        break;
                    }
                }
                
                while ((line = reader.readLine()) != null) {
                    String[] cs = line.split("\t");
                    sequenceSet.add(cs[seqId]);
                }
            } catch (IOException e) {
                LOGGER.error("Error reading sequences for taxonomy analysis", e);
                return false;
            }

            String[] sequences = sequenceSet.toArray(new String[sequenceSet.size()]);
            Arrays.sort(sequences);

            HashMap<String, UnipLCAResult> lcaResultMap = UnipeptRequester.pept2lca(sequences);

            MetaPeptideXMLWriter writer = new MetaPeptideXMLWriter(taxResultFile, MetaConstants.pFind,
                    MetaConstants.labelFree, expNames);

            TaxonomyDatabase td = new TaxonomyDatabase(mrp.getToolPath(ToolType.FUNCDB));

            for (int i = 0; i < sequences.length; i++) {
                if (lcaResultMap.containsKey(sequences[i])) {
                    UnipLCAResult lcaResult = lcaResultMap.get(sequences[i]);
                    int taxId = lcaResult.getTaxon().getId();
                    Taxon lca = td.getTaxonFromId(taxId);

                    if (lca != null) {
                        // Create a simple MetaPeptide for taxonomy analysis
                        double[] intensity = pepIntensityMap.getOrDefault(sequences[i], new double[expNames.length]);
                        PFindPeptide op = new PFindPeptide(sequences[i], 0, 0.0, "", new String[]{sequences[i]}, 
                                new double[]{0.0}, new int[]{1}, new double[][]{intensity});
                        op.setLca(lca);

                        if (lca.getMainParentIds() == null) {
                            lca.setMainParentIds(td.getMainParentTaxonIds(taxId));
                        }

                        int[] mainTaxonIds = lca.getMainParentIds();
                        for (int parentId : mainTaxonIds) {
                            Taxon parentTaxon = td.getTaxonFromId(parentId);
                            if (parentTaxon != null) {
                                if (parentTaxon.getMainParentIds() == null) {
                                    parentTaxon.setMainParentIds(td.getMainParentTaxonIds(parentId));
                                }
                                writer.addTaxon(parentTaxon);
                            }
                        }

                        writer.addPeptide(op);
                    }
                }
            }
            writer.close();

            // Generate taxonomy analysis files (same as DDAIgcPFindTask)
            if (taxResultFile.exists() && taxResultFile.length() > 0) {
                MetaPeptideXMLReader reader = new MetaPeptideXMLReader(taxResultFile);
                File refinedTaxon = new File(taxFile, MetaAlgorithm.Builtin.getName() + ".taxa.refine.csv");
                reader.exportCsv(taxFile, MetaAlgorithm.Unipept, mptp.getLeastPepCount());

                MetaPeptide[] peps = reader.getPeptides();
                Taxon[] taxons = reader.getTaxons();

                File megan = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".biom");
                if (!megan.exists() || megan.length() == 0) {
                    MetaBiomJsonHandler.export(peps, taxons, expNames, megan.getAbsolutePath());
                }

                File allPepTaxa = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".allPepTaxa.csv");
                if (!allPepTaxa.exists() || allPepTaxa.length() == 0) {
                    reader.exportPeptideTaxaAll(allPepTaxa);
                }

                File iMetaLab = new File(taxFile, MetaAlgorithm.Unipept.getName() + ".iMetaLab.tree.csv");
                if (!iMetaLab.exists() || iMetaLab.length() == 0) {
                    MetaTreeHandler.export(peps, taxons, expNames, iMetaLab);
                }

                File jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree.js");
                if (!jsOutput.exists() || jsOutput.length() == 0) {
                    MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
                } else {
                    jsOutput = new File(metaPar.getReportDir() + "\\reports\\metamap\\data\\tree2.js");
                    MetaTreeHandler.exportJS(peps, taxons, expNames, jsOutput);
                }
            }
        }

        updateProgress(95, 100);

        // Functional annotation (same as DDAIgcPFindTask)
        LOGGER.info(getTaskName() + ": functional annotation started");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

        File funcFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
        if (!funcFile.exists()) {
            funcFile.mkdir();
        }

        File final_pro_xml_file = new File(funcFile, "functions.xml");
        if (!final_pro_xml_file.exists() || final_pro_xml_file.length() == 0) {
            // Read protein data for functional annotation
            try (BufferedReader quanProReader = new BufferedReader(new FileReader(final_pro_txt))) {
                String line = quanProReader.readLine(); // Skip header
                // Process protein data for functional annotation
                // This would be similar to DDAIgcPFindTask but using Sage protein data
            } catch (IOException e) {
                LOGGER.error("Error reading protein data for functional annotation", e);
            }
        }

        updateProgress(98, 100);

        // Generate reports
        boolean pepReportOk = false;
        boolean proReportOk = false;
        
        if (final_pep_txt.exists() && final_pep_txt.length() > 0) {
            File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
            if (!report_peptides_summary.exists() || report_peptides_summary.length() == 0) {
                if (summaryMetaFile.exists()) {
                    pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary, summaryMetaFile);
                } else {
                    pepReportOk = task.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt, report_peptides_summary);
                }
            }
        }

        if (final_pro_txt.exists() && final_pro_txt.length() > 0) {
            File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
            if (!report_proteinGroups_summary.exists() || report_proteinGroups_summary.length() == 0) {
                if (summaryMetaFile.exists()) {
                    proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary, summaryMetaFile);
                } else {
                    proReportOk = task.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt, report_proteinGroups_summary);
                }
            }
        }

        // Taxonomy report
        File refinedTaxon = new File(taxFile, MetaAlgorithm.Builtin.getName() + ".taxa.refine.csv");
        File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
        if (!report_taxonomy_summary.exists() || report_taxonomy_summary.length() == 0) {
            if (summaryMetaFile.exists()) {
                task.addTask(MetaPythonReportTask.TAXONOMY, refinedTaxon, report_taxonomy_summary, summaryMetaFile);
            } else {
                task.addTask(MetaPythonReportTask.TAXONOMY, refinedTaxon, report_taxonomy_summary);
            }
        }

        // Functional report
        File funReportTsv = new File(funcFile, "functions_report.tsv");
        File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
        if (!report_function_summary.exists() || report_function_summary.length() == 0) {
            if (summaryMetaFile.exists()) {
                task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary, summaryMetaFile);
            } else {
                task.addTask(MetaPythonReportTask.FUNCTION, funReportTsv, report_function_summary);
            }
        }

        if (!summaryFinished || !pepReportOk || !proReportOk) {
            LOGGER.info(getTaskName() + ": error in generating the reports");
            System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": error in generating the reports");
            return false;
        }

        task.run();

        try {
            if (task.get()) {
                LOGGER.info(getTaskName() + ": report generation successes");
                System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": report generation successes");
                return true;
            }
        } catch (InterruptedException | ExecutionException e) {
            LOGGER.error("Error in report generation", e);
        }

        return false;
    }

    /**
     * Updates protein inference using peptide and PSM data with two-stage filtering.
     * 
     * This method performs the following steps:
     * 1. Builds peptide-to-protein mapping from MetaPSM data
     * 2. Applies protein-level Occam's razor filtering
     * 3. Gets taxonomic information using UnipeptRequester.protinfo()
     * 4. Applies species-level Occam's razor filtering
     * 5. Recalculates protein intensities from peptide data
     * 6. Applies FDR filtering at protein level (< 1%)
     * 7. Updates MetaPeptide.proteins field with final protein names
     * 8. Returns filtered MetaProtein array with FDR < 1%
     * 
     * @param peptides Array of identified peptides
     * @param psms Array of PSMs for protein inference
     * @return Array of filtered MetaProtein objects with FDR < 1%
     */
	protected boolean processPeptideListAndGenerateReports(MetaPeptide[] quanPeps, MetaPSM[] psms,
			HashMap<String, Integer> fileIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName) {

		LOGGER.info(taskName + ": Starting peptide processing and report generation for IGC workflow");

		if (quanPeps == null || quanPeps.length == 0) {
			LOGGER.warn(taskName + ": No quantified peptides provided");
			return false;
		}
		if (psms == null || psms.length == 0) {
			LOGGER.warn(taskName + ": No PSM information provided");
			return false;
		}
		if (fileIdMap == null || fileIdMap.isEmpty()) {
			LOGGER.warn(taskName + ": File mapping is empty");
			return false;
		}

		// First, update proteins to get final filtered protein set
		HashMap<String, UnipProInfo> finalProteinInfoMap = updateProtein(quanPeps);

		if (finalProteinInfoMap.isEmpty()) {
			LOGGER.warn(taskName + ": No proteins survived protein/species refinement");
		}

		// Update peptide protein mappings with final protein names
		String[] razorProteins = new String[quanPeps.length];
		for (int i = 0; i < quanPeps.length; i++) {
			MetaPeptide peptide = quanPeps[i];
			if (peptide == null) {
				continue;
			}
			String[] proteins = peptide.getProteins();
			if (proteins == null) {
				peptide.setProteins(new String[0]);
				continue;
			}
			ArrayList<String> filtered = new ArrayList<String>();
			for (String protein : proteins) {
				if (protein != null && finalProteinInfoMap.containsKey(protein)) {
					filtered.add(protein);
				}
			}
			if (filtered.isEmpty()) {
				peptide.setProteins(new String[0]);
				continue;
			}
			peptide.setProteins(filtered.toArray(new String[filtered.size()]));
			razorProteins[i] = filtered.get(0);
		}

		// Now build PSM maps with the final, filtered protein information
		HashMap<String, MetaPSM> psmMap = new HashMap<String, MetaPSM>();
		HashMap<String, int[]> pepPsmCountMap = new HashMap<String, int[]>();
		HashMap<String, int[]> pepChargeMap = new HashMap<String, int[]>();
		HashMap<String, double[]> pepExpectMap = new HashMap<String, double[]>();
		HashMap<String, Integer> proteinPsmCountMap = new HashMap<String, Integer>();

		int fileCount = fileIdMap.size();
		int[] filePsmCount = new int[fileCount];
		HashSet<String>[] fileSequenceSets = new HashSet[fileCount];
		for (int i = 0; i < fileSequenceSets.length; i++) {
			fileSequenceSets[i] = new HashSet<String>();
		}

		for (MetaPSM psm : psms) {
			if (psm == null) {
				continue;
			}

			String modseq = psm.getModSeq();
			String baseSequence = psm.getBaseSequence();
			String fileName = psm.getFileName();
			Integer fileIdx = fileIdMap.get(fileName);
			if (fileIdx == null) {
				LOGGER.warn(taskName + ": Skip PSM because file '" + fileName + "' is not recognised in file map");
				continue;
			}

			int[] psmCount;
			if (pepPsmCountMap.containsKey(modseq)) {
				psmCount = pepPsmCountMap.get(modseq);
			} else {
				psmCount = new int[fileCount];
				pepPsmCountMap.put(modseq, psmCount);
			}
			psmCount[fileIdx]++;

			filePsmCount[fileIdx]++;
			fileSequenceSets[fileIdx].add(baseSequence);

			int[] charges;
			if (pepChargeMap.containsKey(modseq)) {
				charges = pepChargeMap.get(modseq);
			} else {
				charges = new int[7];
				pepChargeMap.put(modseq, charges);
			}
			int charge = psm.getCharge();
			if (charge >= 0 && charge < charges.length) {
				charges[charge] = 1;
			}

			double[] expect;
			if (pepExpectMap.containsKey(modseq)) {
				expect = pepExpectMap.get(modseq);
						} else {
				expect = new double[fileCount];
				Arrays.fill(expect, 1.0);
				pepExpectMap.put(modseq, expect);
			}
			double pep = psm.getPEP();
			if (pep < expect[fileIdx]) {
				expect[fileIdx] = pep;
			}

			MetaPSM best = psmMap.get(modseq);
			if (best == null || pep < best.getPEP()) {
				psmMap.put(modseq, psm);
			}

			String[] proteins = psm.getProteins();
			if (proteins != null) {
				for (String protein : proteins) {
					if (protein != null && protein.length() > 0) {
						if (proteinPsmCountMap.containsKey(protein)) {
							proteinPsmCountMap.put(protein, proteinPsmCountMap.get(protein) + 1);
						} else {
							proteinPsmCountMap.put(protein, 1);
						}
					}
				}
			}
		}

		for (MetaPeptide peptide : quanPeps) {
			if (peptide == null) {
				continue;
			}
			String modSeq = peptide.getModSeq();
			double[] expect = pepExpectMap.get(modSeq);
			if (expect == null) {
				expect = new double[fileCount];
				Arrays.fill(expect, 1.0);
			}
			peptide.setIndividualPEP(expect);
		}

		// Use the protein info map directly since it already contains the filtered proteins
		HashMap<String, UnipProInfo> proteinInfoMap = finalProteinInfoMap;

		HashMap<String, String> proteinToSpeciesKey = new HashMap<String, String>();
		for (String protein : finalProteinInfoMap.keySet()) {
			UnipProInfo info = proteinInfoMap.get(protein);
			if (info != null) {
				String speciesKey = buildSpeciesKey(info.getTaxId());
				proteinToSpeciesKey.put(protein, speciesKey);
			}
		}

		HashMap<String, double[]> proteinIntensityMap = new HashMap<String, double[]>();
		HashMap<String, Double> proteinTotalIntensityMap = new HashMap<String, Double>();
		HashMap<String, double[]> speciesIntensityMap = new HashMap<String, double[]>();
		HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<String, HashSet<String>>();
		ArrayList<MetaPeptide> peptideList = new ArrayList<MetaPeptide>();

		for (int i = 0; i < quanPeps.length; i++) {
			MetaPeptide peptide = quanPeps[i];
			if (peptide == null) {
				continue;
			}

			String[] proteins = peptide.getProteins();
			if (proteins == null || proteins.length == 0) {
				continue;
			}

			String sequence = peptide.getSequence();
			if (sequence == null) {
				continue;
			}

			int contributingProteins = 0;
			for (String protein : proteins) {
				if (protein != null && finalProteinInfoMap.containsKey(protein)) {
					contributingProteins++;
				}
			}

			if (contributingProteins == 0) {
				continue;
			}

			peptideList.add(peptide);

			if (razorProteins[i] != null) {
				HashSet<String> set;
				if (razorPepCountMap.containsKey(razorProteins[i])) {
					set = razorPepCountMap.get(razorProteins[i]);
			} else {
					set = new HashSet<String>();
					razorPepCountMap.put(razorProteins[i], set);
				}
				set.add(sequence);

				String speciesKey = proteinToSpeciesKey.get(razorProteins[i]);
				if (speciesKey != null) {
					HashSet<String> speciesSet;
					if (razorPepCountMap.containsKey(speciesKey)) {
						speciesSet = razorPepCountMap.get(speciesKey);
					} else {
						speciesSet = new HashSet<String>();
						razorPepCountMap.put(speciesKey, speciesSet);
					}
					speciesSet.add(sequence);
				}
			}

			double[] pepIntensity = new double[peptide.getIntensity().length];
			for (int k = 0; k < pepIntensity.length; k++) {
				pepIntensity[k] = peptide.getIntensity()[k] / (double) contributingProteins;
			}

			for (String protein : proteins) {
				if (protein == null || !finalProteinInfoMap.containsKey(protein)) {
					continue;
				}

				double[] proteinIntensity;
				if (proteinIntensityMap.containsKey(protein)) {
					proteinIntensity = proteinIntensityMap.get(protein);
				} else {
					proteinIntensity = new double[pepIntensity.length];
					proteinIntensityMap.put(protein, proteinIntensity);
				}
				double total = proteinTotalIntensityMap.containsKey(protein)
						? proteinTotalIntensityMap.get(protein)
						: 0.0;

				for (int k = 0; k < pepIntensity.length; k++) {
					proteinIntensity[k] += pepIntensity[k];
					total += pepIntensity[k];
				}
				proteinTotalIntensityMap.put(protein, total);

				HashSet<String> totalSet;
				if (totalPepCountMap.containsKey(protein)) {
					totalSet = totalPepCountMap.get(protein);
				} else {
					totalSet = new HashSet<String>();
					totalPepCountMap.put(protein, totalSet);
				}
				totalSet.add(sequence);

				String speciesKey = proteinToSpeciesKey.get(protein);
				if (speciesKey != null) {
					double[] speciesIntensity;
					if (speciesIntensityMap.containsKey(speciesKey)) {
						speciesIntensity = speciesIntensityMap.get(speciesKey);
			} else {
						speciesIntensity = new double[pepIntensity.length];
						speciesIntensityMap.put(speciesKey, speciesIntensity);
					}
					for (int k = 0; k < pepIntensity.length; k++) {
						speciesIntensity[k] += pepIntensity[k];
					}

					HashSet<String> speciesTotal;
					if (totalPepCountMap.containsKey(speciesKey)) {
						speciesTotal = totalPepCountMap.get(speciesKey);
			} else {
						speciesTotal = new HashSet<String>();
						totalPepCountMap.put(speciesKey, speciesTotal);
					}
					speciesTotal.add(sequence);
				}
			}
		}

		HashMap<String, Integer> filteredProteinPsmCountMap = new HashMap<String, Integer>();
		for (String protein : finalProteinInfoMap.keySet()) {
			if (proteinPsmCountMap.containsKey(protein)) {
				filteredProteinPsmCountMap.put(protein, proteinPsmCountMap.get(protein));
			}
		}
		proteinPsmCountMap = filteredProteinPsmCountMap;

		HashMap<String, Integer> speciesPsmCountMap = new HashMap<String, Integer>();
		for (Map.Entry<String, Integer> entry : proteinPsmCountMap.entrySet()) {
			String speciesKey = proteinToSpeciesKey.get(entry.getKey());
			if (speciesKey != null) {
				if (speciesPsmCountMap.containsKey(speciesKey)) {
					speciesPsmCountMap.put(speciesKey, speciesPsmCountMap.get(speciesKey) + entry.getValue());
			} else {
					speciesPsmCountMap.put(speciesKey, entry.getValue());
				}
			}
		}

		String[] proteins = proteinTotalIntensityMap.keySet().toArray(new String[proteinTotalIntensityMap.size()]);
		Arrays.sort(proteins, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				double inten1 = proteinTotalIntensityMap.get(o1);
				double inten2 = proteinTotalIntensityMap.get(o2);
				if (inten1 > inten2)
					return -1;
				if (inten1 < inten2)
					return 1;
				return 0;
			}
		});

		HashMap<String, Integer> proteinIdMap = new HashMap<String, Integer>();
		for (int i = 0; i < proteins.length; i++) {
			proteinIdMap.put(proteins[i], i);
		}

		if (!generateIdentificationSummary(fileIdMap, filePsmCount, fileSequenceSets, taskName)) {
			return false;
		}

		updateProgress(85, 100);

		if (!generatePeptideReport(peptideList, psmMap, pepPsmCountMap, pepChargeMap, proteinIdMap, fileList, quanMode,
				intensityTitle, taskName)) {
			return false;
		}

		if (!generateProteinReport(proteins, proteinIntensityMap, proteinTotalIntensityMap, totalPepCountMap,
					razorPepCountMap, proteinPsmCountMap, proteinIdMap, quanMode, intensityTitle, taskName)) {
			return false;
		}

		updateProgress(89, 100);

		if (!performFunctionalAnnotation(proteins, proteinIntensityMap, proteinTotalIntensityMap, totalPepCountMap, proteinInfoMap, intensityTitle, taskName)) {
			return false;
		}

		if (!performTaxonomyAnalysis(proteins, proteinIntensityMap, proteinPsmCountMap, totalPepCountMap, razorPepCountMap,
				proteinInfoMap, intensityTitle, taskName)) {
			return false;
		}

		updateProgress(95, 100);

		if (!generateFinalReports(taskName)) {
			return false;
		}

		LOGGER.info(taskName + ": Peptide processing and report generation completed successfully for IGC workflow");
		return true;
	}
	
	protected HashMap<String, UnipProInfo> updateProtein(MetaPeptide[] peptides) {
		LOGGER.info(getTaskName() + ": Starting protein inference and filtering using peptide and PSM data");

		// Step 1: Build peptide-to-protein mapping from MetaPeptide data
		HashMap<String, HashSet<String>> peptideToProteins = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> proteinToPeptides = new HashMap<String, HashSet<String>>();

		for (MetaPeptide peptide : peptides) {
			if (peptide != null && peptide.getSequence() != null && peptide.getProteins() != null) {
				String sequence = peptide.getSequence();
				String[] peptideProteins = peptide.getProteins();

				if (!peptideToProteins.containsKey(sequence)) {
					peptideToProteins.put(sequence, new HashSet<String>());
				}

				for (String protein : peptideProteins) {
					if (protein != null && !protein.isEmpty()) {
						peptideToProteins.get(sequence).add(protein);

						if (!proteinToPeptides.containsKey(protein)) {
							proteinToPeptides.put(protein, new HashSet<String>());
						}
						proteinToPeptides.get(protein).add(sequence);
					}
				}
			}
		}

		LOGGER.info(getTaskName() + ": Found " + proteinToPeptides.size() + " proteins with " + peptideToProteins.size()
				+ " unique peptides");

		// Step 2: Extract unique peptide sequences for Unipept API call
		HashSet<String> uniqueSequences = new HashSet<String>(peptideToProteins.keySet());
		if (uniqueSequences.isEmpty()) {
			LOGGER.warn(getTaskName() + ": No valid peptide sequences found for protein inference");
			return new HashMap<String, UnipProInfo>();
		}

		String[] sequences = uniqueSequences.toArray(new String[uniqueSequences.size()]);
		LOGGER.info(getTaskName() + ": Calling Unipept API for " + sequences.length + " unique peptide sequences");

		// Step 3: Call UnipeptRequester.pept2prot to get protein mappings
		HashMap<String, UnipProResult> pept2protMap = UnipeptRequester.pept2prot(sequences);

		if (pept2protMap == null || pept2protMap.isEmpty()) {
			LOGGER.warn(getTaskName() + ": No protein mappings returned from Unipept API");
			return new HashMap<String, UnipProInfo>();
		}

		LOGGER.info(getTaskName() + ": Received " + pept2protMap.size() + " protein mappings from Unipept API");

		// Step 4: Apply protein-level Occam's razor (similar to DDA MAG workflow)
		HashSet<String> selectedProteins = applyOccamsRazor(proteinToPeptides, uniqueSequences, 0.95);
		LOGGER.info(getTaskName() + ": Protein-level filtering selected " + selectedProteins.size() + " proteins");

		// Step 5: Get taxonomic information for selected proteins using protinfo
		HashSet<String> uniprotProteins = new HashSet<String>();
		for (String peptide : selectedProteins) {
			HashSet<String> peptideProteins = proteinToPeptides.get(peptide);
			if (peptideProteins != null) {
				uniprotProteins.addAll(peptideProteins);
			}
		}

		// Get protein information from Unipept
		String[] proteinArray = uniprotProteins.toArray(new String[uniprotProteins.size()]);
		HashMap<String, UnipProInfo> proteinInfoMap = UnipeptRequester.protinfo(proteinArray);

		if (proteinInfoMap == null || proteinInfoMap.isEmpty()) {
			LOGGER.warn(getTaskName() + ": No protein information returned from Unipept protinfo API");
			return new HashMap<String, UnipProInfo>();
		}

		// Step 6: Build protein-to-species mapping
		HashMap<String, Integer> proteinToSpeciesMap = new HashMap<String, Integer>();
		HashMap<Integer, HashSet<String>> speciesToProteins = new HashMap<Integer, HashSet<String>>();

		for (String protein : selectedProteins) {
			UnipProInfo info = proteinInfoMap.get(protein);
			if (info != null) {
				int speciesId = info.getTaxId(); // Already at species level
				proteinToSpeciesMap.put(protein, speciesId);

				if (!speciesToProteins.containsKey(speciesId)) {
					speciesToProteins.put(speciesId, new HashSet<String>());
				}
				speciesToProteins.get(speciesId).add(protein);
			}
		}

		LOGGER.info(getTaskName() + ": Mapped " + proteinToSpeciesMap.size() + " proteins to "
				+ speciesToProteins.size() + " species");

		// Step 7: Apply species-level Occam's razor based on peptide coverage
		// Build species-to-peptides mapping for Occam's razor
		HashMap<String, HashSet<String>> speciesToPeptidesString = new HashMap<String, HashSet<String>>();
		HashMap<String, Integer> stringToIntegerSpecies = new HashMap<String, Integer>();
		
		for (Map.Entry<Integer, HashSet<String>> entry : speciesToProteins.entrySet()) {
			Integer speciesId = entry.getKey();
			HashSet<String> speciesProteins = entry.getValue();
			String speciesKey = String.valueOf(speciesId);
			stringToIntegerSpecies.put(speciesKey, speciesId);
			
			// Build peptide set for this species
			HashSet<String> speciesPeptides = new HashSet<String>();
			for (String protein : speciesProteins) {
				HashSet<String> proteinPeptides = proteinToPeptides.get(protein);
				if (proteinPeptides != null) {
					speciesPeptides.addAll(proteinPeptides);
				}
			}
			speciesToPeptidesString.put(speciesKey, speciesPeptides);
		}
		
		// Use total unique peptides for threshold calculation
		HashSet<String> selectedSpeciesString = applyOccamsRazor(speciesToPeptidesString, uniqueSequences, 0.95);
		
		// Convert back to Integer set
		HashSet<Integer> selectedSpecies = new HashSet<Integer>();
		for (String speciesKey : selectedSpeciesString) {
			selectedSpecies.add(stringToIntegerSpecies.get(speciesKey));
		}
		LOGGER.info(getTaskName() + ": Species-level filtering selected " + selectedSpecies.size() + " species");

		// Step 8: Update protein list to remove proteins from discarded species
		HashSet<String> finalSelectedProteins = new HashSet<String>();
		for (String protein : selectedProteins) {
			Integer speciesId = proteinToSpeciesMap.get(protein);
			if (speciesId != null && selectedSpecies.contains(speciesId)) {
				finalSelectedProteins.add(protein);
			}
		}

		LOGGER.info(getTaskName() + ": Final protein set contains " + finalSelectedProteins.size()
				+ " proteins after species filtering");

		// Step 9: Apply FDR filtering to protein names (no MetaProtein objects created here)
		// Note: Intensity calculations are done in processPeptideListAndGenerateReports method
		ArrayList<String> proteinList = new ArrayList<String>(finalSelectedProteins);

		// Sort proteins by peptide count (descending) for FDR calculation
		proteinList.sort((p1, p2) -> {
			int count1 = proteinToPeptides.get(p1) != null ? proteinToPeptides.get(p1).size() : 0;
			int count2 = proteinToPeptides.get(p2) != null ? proteinToPeptides.get(p2).size() : 0;
			return Integer.compare(count2, count1);
		});

		int decoyCount = 0;
		int targetCount = 0;
		HashSet<String> filteredProteins = new HashSet<String>();

		for (String protein : proteinList) {
			boolean isDecoy = protein.startsWith("REV_") || protein.contains("DECOY") || protein.startsWith("REV")
					|| protein.endsWith("_REV");

			if (isDecoy) {
				decoyCount++;
			} else {
				targetCount++;
			}

			// Calculate FDR
			double fdr = (double) decoyCount / (double) (targetCount + decoyCount);

			if (fdr < 0.01) { // FDR < 1%
				filteredProteins.add(protein);
					} else {
				// Stop when FDR >= 1%
				break;
			}
		}

		LOGGER.info(
				getTaskName() + ": FDR filtering complete - " + filteredProteins.size() + " proteins kept at FDR < 1%");

		// Step 11: Update peptide protein mappings with final protein names
		updatePeptideProteinMappings(peptides, filteredProteins);

		// Step 12: Build final protein info map with only filtered proteins
		HashMap<String, UnipProInfo> finalProteinInfoMap = new HashMap<String, UnipProInfo>();
		for (String protein : filteredProteins) {
			UnipProInfo info = proteinInfoMap.get(protein);
			if (info != null) {
				finalProteinInfoMap.put(protein, info);
			}
		}

		LOGGER.info(getTaskName() + ": Protein inference complete - returning " + finalProteinInfoMap.size()
				+ " filtered proteins with taxonomic information");

		return finalProteinInfoMap;
	}

	/**
	 * Gets detailed protein information for UniProt proteins using the protinfo API.
	 * This method extracts UniProt protein IDs from MetaProtein objects and retrieves
	 * detailed information including protein names, taxonomic data, and functional annotations.
	 * 
	 * @param proteins Array of MetaProtein objects with UniProt protein names
	 * @return HashMap mapping UniProt protein IDs to their detailed information
	 */
	protected HashMap<String, UnipProInfo> getProteinInfo(MetaProtein[] proteins) {
		LOGGER.info(getTaskName() + ": Starting protein information retrieval using Unipept protinfo API");
		
		if (proteins == null || proteins.length == 0) {
			LOGGER.warn(getTaskName() + ": No proteins provided for information retrieval");
			return new HashMap<String, UnipProInfo>();
		}
		
		// Step 1: Extract unique UniProt protein IDs from MetaProtein array
		HashSet<String> uniqueProteinIds = new HashSet<String>();
		for (MetaProtein protein : proteins) {
			if (protein != null && protein.getName() != null && !protein.getName().isEmpty()) {
				// Only include non-decoy proteins (UniProt IDs)
				boolean isDecoy = protein.getName().startsWith("REV_") || protein.getName().contains("DECOY") ||
								protein.getName().startsWith("REV") || protein.getName().endsWith("_REV");
				
				if (!isDecoy) {
					uniqueProteinIds.add(protein.getName());
				}
			}
		}
		
		if (uniqueProteinIds.isEmpty()) {
			LOGGER.warn(getTaskName() + ": No valid UniProt protein IDs found for information retrieval");
			return new HashMap<String, UnipProInfo>();
		}
		
		// Convert HashSet to String array for Unipept API call
		String[] proteinIds = uniqueProteinIds.toArray(new String[uniqueProteinIds.size()]);
		LOGGER.info(getTaskName() + ": Calling Unipept protinfo API for " + proteinIds.length + " unique protein IDs");
		
		// Step 2: Call UnipeptRequester.protinfo to get protein information
		HashMap<String, UnipProInfo> proteinInfoMap = UnipeptRequester.protinfo(proteinIds);
		
		if (proteinInfoMap == null || proteinInfoMap.isEmpty()) {
			LOGGER.warn(getTaskName() + ": No protein information returned from Unipept protinfo API");
			return new HashMap<String, UnipProInfo>();
		}
		
		LOGGER.info(getTaskName() + ": Successfully retrieved information for " + proteinInfoMap.size() + 
				   " proteins from Unipept protinfo API");
		
		// Log some sample protein information for verification
		int sampleCount = 0;
		for (String proteinId : proteinInfoMap.keySet()) {
			if (sampleCount < 3) { // Log first 3 proteins as examples
				UnipProInfo info = proteinInfoMap.get(proteinId);
				LOGGER.info(getTaskName() + ": Sample protein info - ID: " + proteinId + 
						   ", Name: " + info.getName() + ", Taxon ID: " + info.getTaxId() +
						   ", EC numbers: " + java.util.Arrays.toString(info.getEc()) +
						   ", GO terms: " + java.util.Arrays.toString(info.getGo()));
				sampleCount++;
			}
		}
		
		return proteinInfoMap;
	}


	/**
	 * Exports functional annotation results for proteins.
	 * This method handles the common functional annotation workflow for IGC tasks.
	 */
    protected MetaProteinAnnoIgc[] funcAnnoSqlite(MetaProtein[] pros, String[] fileNames) {
        // This method would contain the common functional annotation logic
        // Implementation would be similar to MetaPilotMagTask but adapted for IGC workflows
        
        LOGGER.info(getTaskName() + ": functional annotation started");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": functional annotation started");

        // TODO: Implement functional annotation logic specific to IGC workflows
        // This would include COG, NOG, KEGG, GO, and EC annotations
        
        return new MetaProteinAnnoIgc[0];
    }

    /**
     * Exports summary statistics for the analysis.
     * This method generates the ID summary report common to all IGC workflows.
     */
    protected boolean exportSummary(MetaPythonReportTask task) {

        File final_summary_txt = new File(this.metaPar.getDbSearchResultFile(), "final_summary.tsv");
        try {

            PrintWriter countWriter = new PrintWriter(final_summary_txt);
            StringBuilder titlesb = new StringBuilder();
            titlesb.append("Raw file").append("\t");
            titlesb.append("Experiment").append("\t");
            titlesb.append("MS/MS").append("\t");
            titlesb.append("MS/MS Identified").append("\t");
            titlesb.append("MS/MS Identified [%]").append("\t");
            titlesb.append("Peptide Sequences Identified");

            countWriter.println(titlesb);

            int totalms2count = 0;
            int totalpsmCount = 0;
            int totalUnipepCount = 0;

            BufferedReader idReader = new BufferedReader(
                    new FileReader(new File(metaPar.getDbSearchResultFile(), "id.tsv")));
            HashMap<String, int[]> idMap = new HashMap<String, int[]>();
            String line = idReader.readLine();
            while ((line = idReader.readLine()) != null) {
                String[] cs = line.split("\t");
                idMap.put(cs[0], new int[] { Integer.parseInt(cs[1]), Integer.parseInt(cs[2]) });
            }
            idReader.close();

            File ms2CountFile = new File(spectraDirFile, "ms2Count.txt");
            BufferedReader ms2CountReader = new BufferedReader(new FileReader(ms2CountFile));
            line = ms2CountReader.readLine();

            while ((line = ms2CountReader.readLine()) != null) {
                String[] cs = line.split("\t");
                String rawName = cs[0];

                if (idMap.containsKey(rawName)) {

                    StringBuilder sb = new StringBuilder();
                    sb.append(cs[0]).append("\t");

                    if (this.expNameMap.containsKey(cs[0])) {
                        sb.append(expNameMap.get(cs[0])).append("\t");
                    } else {
                        sb.append(cs[0]).append("\t");
                    }

                    int ms2count = Integer.parseInt(cs[1]);

                    int[] idCounts = idMap.get(rawName);

                    sb.append(ms2count).append("\t");
                    sb.append(idCounts[1]).append("\t");

                    double ratio = 0;
                    if (ms2count > 0) {
                        ratio = (double) idCounts[1] / (double) ms2count * 100.0;
                    }

                    sb.append(FormatTool.getDF2().format(ratio)).append("\t");

                    sb.append(idCounts[0]);

                    countWriter.println(sb);

                    totalms2count += ms2count;
                    totalpsmCount += idCounts[1];
                    totalUnipepCount += idCounts[0];
                }
            }

            ms2CountReader.close();

            StringBuilder sb = new StringBuilder();
            sb.append("Total").append("\t");
            sb.append("\t");
            sb.append(totalms2count).append("\t");
            sb.append(totalpsmCount).append("\t");

            double ratio = 0;
            if (totalms2count > 0) {
                ratio = (double) totalpsmCount / (double) totalms2count * 100.0;
            }

            sb.append(FormatTool.getDF2().format(ratio)).append("\t");

            sb.append(totalUnipepCount);

            countWriter.print(sb);

            countWriter.close();

        } catch (IOException e) {
            // TODO Auto-generated catch block
            LOGGER.error(getTaskName() + ": error in writing PSM information to " + final_summary_txt.getName(), e);
            System.err.println(format.format(new Date()) + "\t" + getTaskName()
                    + ": error in writing PSM information to " + final_summary_txt.getName());
        }

        if (!final_summary_txt.exists() || final_summary_txt.length() == 0) {
            return false;

        } else {
            boolean finish = false;
            File report_ID_summary = new File(metaPar.getReportDir(), "report_ID_summary.html");
            if (summaryMetaFile != null && summaryMetaFile.exists()) {
                finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary,
                        summaryMetaFile);
            } else {
                finish = task.addTask(MetaPythonReportTask.ID_SUMMARY, final_summary_txt, report_ID_summary);
            }

            return finish;
        }
    }

    /**
     * Helper method to refine proteins using Occam's razor principle (similar to DDA MAG workflow).
     * 
     * @param proteinToPeptides Map from protein names to their associated peptides
     * @param totalPeptides Total number of unique peptides
     * @param threshold Threshold for protein selection (e.g., 0.95 for 95% peptide coverage)
     * @return Set of selected protein names
     */
	private boolean generateIdentificationSummary(HashMap<String, Integer> fileIdMap, int[] filePsmCount,
			HashSet<String>[] fileSequenceSets, String taskName) {

		try (PrintWriter idWriter = new PrintWriter(new File(metaPar.getDbSearchResultFile(), "id.tsv"))) {
			idWriter.println("Raw file\tUnique peptide count\tPSM count");

			String[] rawNames = fileIdMap.keySet().toArray(new String[fileIdMap.size()]);
			Arrays.sort(rawNames, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return fileIdMap.get(o1) - fileIdMap.get(o2);
				}
			});

			for (int i = 0; i < rawNames.length; i++) {
				idWriter.println(rawNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
			}
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing identification information", e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing identification information");
			return false;
		}
		return true;
	}

	private boolean generatePeptideReport(ArrayList<MetaPeptide> peptideList, HashMap<String, MetaPSM> psmMap,
			HashMap<String, int[]> pepPsmCountMap, HashMap<String, int[]> pepChargeMap,
			HashMap<String, Integer> proteinIdMap, ArrayList<String> fileList, String quanMode, String[] intensityTitle,
			String taskName) {

		LOGGER.info(taskName + ": exporting peptide report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting peptide report started");

		try {
			final_pep_txt = new File(this.metaPar.getDbSearchResultFile(), "final_peptides.tsv");
			PrintWriter writer = new PrintWriter(final_pep_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Sequence").append("\t");
			titlesb.append("Base sequence").append("\t");
			titlesb.append("Length").append("\t");
			titlesb.append("Missed cleavages").append("\t");
			titlesb.append("Mass").append("\t");
			titlesb.append("Proteins").append("\t");
			titlesb.append("Charges").append("\t");
			titlesb.append("Score").append("\t");
			titlesb.append("PEP").append("\t");

			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("Identification type " + fileList.get(i)).append("\t");
			}

			for (int i = 0; i < fileList.size(); i++) {
				titlesb.append("MS/MS count " + fileList.get(i)).append("\t");
			}

			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < intensityTitle.length; i++) {
				titlesb.append(intensityTitle[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");
			titlesb.append("Protein group IDs").append("\t");
			titlesb.append("Total MS/MS count");

			writer.println(titlesb);

			int id = 0;
			for (MetaPeptide peptide : peptideList) {
				String fullSeq = peptide.getModSeq();
				MetaPSM psm = psmMap.get(fullSeq);

				StringBuilder sb = new StringBuilder();
				sb.append(peptide.getSequence()).append("\t");
				sb.append(peptide.getSequence()).append("\t");
				sb.append(peptide.getLength()).append("\t");
				sb.append(peptide.getMissCleave()).append("\t");

				double mass = peptide.getMass();
				if (mass == 0 && psm != null) {
					mass = psm.getCalMass();
				}
				sb.append(mass).append("\t");

				String[] pros = peptide.getProteins();
				if (pros != null) {
					StringBuilder proSb = new StringBuilder();
					for (String pro : pros) {
						proSb.append(pro).append(";");
					}
					if (proSb.length() > 0) {
						proSb.deleteCharAt(proSb.length() - 1);
					}
					sb.append(proSb);
				}
				sb.append("\t");

				int[] charge = pepChargeMap.get(fullSeq);
				if (charge != null) {
					StringBuilder chargesb = new StringBuilder();
					for (int j = 0; j < charge.length; j++) {
						if (charge[j] == 1) {
							chargesb.append(j).append(";");
						}
					}
					if (chargesb.length() > 0) {
						if (chargesb.charAt(chargesb.length() - 1) == ';') {
							chargesb.deleteCharAt(chargesb.length() - 1);
						}
						sb.append(chargesb).append("\t");
					} else {
						sb.append("2\t");
					}
				} else {
					sb.append("2\t");
				}

				sb.append(peptide.getScore()).append("\t");
				sb.append(peptide.getPEP()).append("\t");

				int totalMsCount = 0;
				int[] psmCount = pepPsmCountMap.get(fullSeq);
				if (psmCount != null) {
					double[] intensity = peptide.getIntensity();
					for (int j = 0; j < psmCount.length; j++) {
						if (psmCount[j] > 0) {
							totalMsCount += psmCount[j];
							sb.append("By MS/MS").append("\t");
						} else {
							if (intensity[j] > 0) {
								sb.append("MBR\t");
							} else {
								sb.append("NotDetected\t");
							}
						}
					}

					for (int j = 0; j < psmCount.length; j++) {
						sb.append(psmCount[j]).append("\t");
					}
				} else {
					for (int j = 0; j < fileList.size(); j++) {
						sb.append("NotDetected\t");
					}
					for (int j = 0; j < fileList.size(); j++) {
						sb.append(0).append("\t");
					}
				}

				double totalIntensity = 0;
				double[] intensity = peptide.getIntensity();
				for (int j = 0; j < intensity.length; j++) {
					totalIntensity += intensity[j];
				}

				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				for (int j = 0; j < intensity.length; j++) {
					sb.append(BigInteger.valueOf((long) Math.floor(intensity[j]))).append("\t");
				}

				sb.append("").append("\t");
				sb.append("").append("\t");
				sb.append(id).append("\t");

				StringBuilder groupSb = new StringBuilder();
				for (String pro : pros) {
					if (proteinIdMap.containsKey(pro)) {
						groupSb.append(proteinIdMap.get(pro)).append(";");
					}
				}
				if (groupSb.length() > 0) {
					groupSb.deleteCharAt(groupSb.length() - 1);
				}
				sb.append(groupSb).append("\t");

				sb.append(totalMsCount);

				writer.println(sb);
				id++;
			}

			writer.close();

		} catch (Exception e) {
			LOGGER.error(taskName + ": error in writing the peptides to " + final_pep_txt, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the peptides to " + final_pep_txt);
			return false;
		}

		LOGGER.info(taskName + ": peptide report has been exported to " + final_pep_txt);
		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": peptide report has been exported to " + final_pep_txt);

		if (!final_pep_txt.exists() || final_pep_txt.length() == 0) {
			return false;
		}

		return true;
	}

	private boolean generateProteinReport(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntensityMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> proteinPsmCountMap,
			HashMap<String, Integer> proteinIdMap, String quanMode, String[] intensityTitle, String taskName) {

		LOGGER.info(taskName + ": exporting protein report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting protein report started");

		try {
			final_pro_txt = new File(this.metaPar.getDbSearchResultFile(), "final_proteins.tsv");
			PrintWriter writer = new PrintWriter(final_pro_txt);

			StringBuilder titlesb = new StringBuilder();
			titlesb.append("Protein IDs").append("\t");
			titlesb.append("Peptide counts (all)").append("\t");
			titlesb.append("Peptide counts (razor)").append("\t");
			titlesb.append("Number of PSMs").append("\t");
			titlesb.append("PEP").append("\t");
			titlesb.append("Intensity").append("\t");

			for (int i = 0; i < intensityTitle.length; i++) {
				titlesb.append(intensityTitle[i]).append("\t");
			}

			titlesb.append("Reverse").append("\t");
			titlesb.append("Potential contaminant").append("\t");
			titlesb.append("id").append("\t");

			writer.println(titlesb);

			for (String protein : proteins) {
				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");

				int totalPepCount = 0;
				if (totalPepCountMap.containsKey(protein)) {
					totalPepCount = totalPepCountMap.get(protein).size();
				}
				sb.append(totalPepCount).append("\t");

				int razorPepCount = 0;
				if (razorPepCountMap.containsKey(protein)) {
					razorPepCount = razorPepCountMap.get(protein).size();
				}
				sb.append(razorPepCount).append("\t");

				int psmCount = 0;
				if (proteinPsmCountMap.containsKey(protein)) {
					psmCount = proteinPsmCountMap.get(protein);
				}
				sb.append(psmCount).append("\t");

				double totalIntensity = proteinTotalIntensityMap.containsKey(protein)
						? proteinTotalIntensityMap.get(protein)
						: 0.0;
				sb.append(0.0).append("\t");
				sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

				double[] intensity = proteinIntensityMap.get(protein);
				if (intensity != null) {
					for (int i = 0; i < intensity.length; i++) {
						sb.append(BigInteger.valueOf((long) Math.floor(intensity[i]))).append("\t");
					}
				} else {
					for (int i = 0; i < intensityTitle.length; i++) {
						sb.append(0).append("\t");
					}
				}

				sb.append("").append("\t");
				sb.append("").append("\t");
				sb.append(proteinIdMap.get(protein));

				writer.println(sb);
			}

			writer.close();

		} catch (Exception e) {
			LOGGER.error(taskName + ": error in writing the proteins to " + final_pro_txt, e);
			System.err.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the proteins to " + final_pro_txt);
			return false;
		}

		LOGGER.info(taskName + ": protein report has been exported to " + final_pro_txt);
		System.out.println(
				format.format(new Date()) + "\t" + taskName + ": protein report has been exported to " + final_pro_txt);

		if (!final_pro_txt.exists() || final_pro_txt.length() == 0) {
			return false;
		}

		return true;
	}

	private boolean performFunctionalAnnotation(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Double> proteinTotalIntensityMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, UnipProInfo> proteinInfoMap, String[] intensityTitle, String taskName) {

		LOGGER.info(taskName + ": functional annotation started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation started");

		File funcFile = new File(metaPar.getDbSearchResultFile(), "functional_annotation");
		if (!funcFile.exists()) {
			funcFile.mkdir();
		}

		// Write functions.tsv directly from UnipProInfo data (no XML generation needed)
		File funTsv = new File(funcFile, "functions.tsv");
		try (PrintWriter writer = new PrintWriter(funTsv)) {
			// Build header
			StringBuilder header = new StringBuilder();
			header.append("Protein ID\t");
			header.append("Protein Name\t");
			header.append("EC Number\t");
			header.append("GO Terms\t");
			header.append("Taxonomic ID\t");
			header.append("Peptide Count\t");
			header.append("Total Intensity");
			
			// Add intensity columns
			for (String intensity : intensityTitle) {
				if (intensity.startsWith("Intensity ")) {
					header.append("\t").append(intensity);
				} else {
					header.append("\tIntensity ").append(intensity);
				}
			}
			
			writer.println(header);

			// Write protein data
			for (String protein : proteins) {
				UnipProInfo info = proteinInfoMap.get(protein);
				if (info == null) {
					continue; // Skip proteins without UnipProInfo
				}

				StringBuilder sb = new StringBuilder();
				sb.append(protein).append("\t");
				sb.append(info.getName() != null ? info.getName() : "").append("\t");
				
				// Extract EC numbers (if available)
				String ecNumbers = "";
				if (info.getEc() != null && info.getEc().length > 0) {
					ecNumbers = String.join(";", info.getEc());
				}
				sb.append(ecNumbers).append("\t");
				
				// Extract GO terms (if available)
				String goTerms = "";
				if (info.getGo() != null && info.getGo().length > 0) {
					goTerms = String.join(";", info.getGo());
				}
				sb.append(goTerms).append("\t");
				
				sb.append(info.getTaxId()).append("\t");
				
				// Add peptide count and intensity data
				if (totalPepCountMap.containsKey(protein)) {
					sb.append(totalPepCountMap.get(protein).size()).append("\t");
				} else {
					sb.append("0").append("\t");
				}
				
				if (proteinTotalIntensityMap.containsKey(protein)) {
					sb.append(proteinTotalIntensityMap.get(protein));
				} else {
					sb.append("0");
				}
				
				// Add individual intensity values
				if (proteinIntensityMap.containsKey(protein)) {
					double[] intensities = proteinIntensityMap.get(protein);
					for (double intensity : intensities) {
						sb.append("\t").append(intensity);
					}
				}
				
				writer.println(sb);
			}

		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing functional annotation to " + funTsv.getName(), e);
			System.err.println(format.format(new Date()) + "\t" + taskName 
					+ ": error in writing functional annotation to " + funTsv.getName());
			return false;
		}

		LOGGER.info(taskName + ": functional annotation completed - wrote " + proteins.length + " proteins to " + funTsv.getName());
		System.out.println(format.format(new Date()) + "\t" + taskName + ": functional annotation completed");

		return true;
	}

	private boolean performTaxonomyAnalysis(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, Integer> proteinPsmCountMap, HashMap<String, HashSet<String>> totalPepCountMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, UnipProInfo> proteinInfoMap,
			String[] intensityTitle, String taskName) {

		LOGGER.info(taskName + ": taxonomy analysis started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis started");

		File taxFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis");
		if (!taxFile.exists()) {
			taxFile.mkdir();
		}

		// Build protein-to-taxon mapping from UnipProInfo
		HashMap<String, Integer> proteinToTaxonMap = new HashMap<String, Integer>();
		HashMap<Integer, HashSet<String>> taxonToProteins = new HashMap<Integer, HashSet<String>>();

		for (String protein : proteins) {
			UnipProInfo info = proteinInfoMap.get(protein);
			if (info != null && info.getTaxId() > 0) {
				int taxonId = info.getTaxId();
				proteinToTaxonMap.put(protein, taxonId);

				if (!taxonToProteins.containsKey(taxonId)) {
					taxonToProteins.put(taxonId, new HashSet<String>());
				}
				taxonToProteins.get(taxonId).add(protein);
			}
		}

		LOGGER.info(taskName + ": Mapped " + proteinToTaxonMap.size() + " proteins to " + taxonToProteins.size()
				+ " taxonomic IDs");

		// Get taxonomic lineage information for each taxon ID (similar to MAG workflow)
		HashMap<Integer, Taxon> taxonInfoMap = new HashMap<Integer, Taxon>();
		TaxonomyDatabaseSQL taxDb = null;
		try {
			taxDb = new TaxonomyDatabaseSQL(mrp.getToolPath(ToolType.FUNCDB));
			for (Integer taxonId : taxonToProteins.keySet()) {
				Taxon taxon = taxDb.getTaxonFromId(taxonId);
				if (taxon != null) {
					taxonInfoMap.put(taxonId, taxon);
				}
			}
		} catch (Exception e) {
			LOGGER.error(taskName + ": Error accessing taxonomy database", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": Error accessing taxonomy database");
			return false;
		}

		// Build taxonomic rank maps (similar to MAG workflow but protein-based)
		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		HashMap<String, String[]>[] taxRankMaps = new HashMap[mainRanks.length];
		HashMap<String, double[]>[] taxIntensityMaps = new HashMap[mainRanks.length];
		HashMap<String, Integer>[] taxPsmCountMaps = new HashMap[mainRanks.length];
		HashMap<String, HashSet<String>>[] taxPepCountMaps = new HashMap[mainRanks.length];
		HashMap<String, Integer>[] taxProteinCountMaps = new HashMap[mainRanks.length];

		for (int i = 0; i < taxRankMaps.length; i++) {
			taxRankMaps[i] = new HashMap<String, String[]>();
			taxIntensityMaps[i] = new HashMap<String, double[]>();
			taxPsmCountMaps[i] = new HashMap<String, Integer>();
			taxPepCountMaps[i] = new HashMap<String, HashSet<String>>();
			taxProteinCountMaps[i] = new HashMap<String, Integer>();
		}

		// Process each protein and aggregate to taxonomic ranks (like MAG workflow)
		for (String protein : proteins) {
			if (!proteinToTaxonMap.containsKey(protein)) {
				continue; // Skip proteins without taxonomic information
			}

			Integer taxonId = proteinToTaxonMap.get(protein);
			Taxon taxon = taxonInfoMap.get(taxonId);
			if (taxon == null) {
				continue;
			}

			// Get protein data
			double[] proteinIntensity = proteinIntensityMap.get(protein);
			Integer psmCount = proteinPsmCountMap.get(protein);
			HashSet<String> pepCount = totalPepCountMap.get(protein);

			if (proteinIntensity == null || psmCount == null || pepCount == null) {
				continue;
			}

			// Get taxonomic lineage using the same approach as existing IGC workflow
			if (taxon.getMainParentIds() == null) {
				try {
					taxon.setMainParentIds(taxDb.getMainParentTaxonIds(taxonId));
				} catch (Exception e) {
					LOGGER.warn(taskName + ": Could not get lineage for taxon " + taxonId);
					continue;
				}
			}

			int[] mainTaxonIds = taxon.getMainParentIds();
			if (mainTaxonIds == null || mainTaxonIds.length == 0) {
				continue;
			}

			// Process each taxonomic rank in the lineage
			for (int rankIndex = 0; rankIndex < mainTaxonIds.length && rankIndex < mainRanks.length; rankIndex++) {
				int parentTaxonId = mainTaxonIds[rankIndex];
				Taxon parentTaxon = taxonInfoMap.get(parentTaxonId);
				if (parentTaxon == null) {
					// Get parent taxon if not already cached (reuse existing taxDb)
					try {
						parentTaxon = taxDb.getTaxonFromId(parentTaxonId);
						if (parentTaxon != null) {
							taxonInfoMap.put(parentTaxonId, parentTaxon);
						}
					} catch (Exception e) {
						continue;
					}
				}

				if (parentTaxon == null) {
					continue;
				}

				String taxonName = parentTaxon.getName();
				if (taxonName == null || taxonName.isEmpty()) {
					continue;
				}

				// Build lineage for this rank
				String[] lineage = new String[rankIndex + 1];
				for (int i = 0; i <= rankIndex && i < mainTaxonIds.length; i++) {
					int lineageTaxonId = mainTaxonIds[i];
					Taxon lineageTaxon = taxonInfoMap.get(lineageTaxonId);
					if (lineageTaxon != null) {
						lineage[i] = lineageTaxon.getName();
					}
				}

				if (!taxRankMaps[rankIndex].containsKey(taxonName)) {
					taxRankMaps[rankIndex].put(taxonName, lineage);
					taxIntensityMaps[rankIndex].put(taxonName, new double[proteinIntensity.length]);
					taxPsmCountMaps[rankIndex].put(taxonName, 0);
					taxPepCountMaps[rankIndex].put(taxonName, new HashSet<String>());
					taxProteinCountMaps[rankIndex].put(taxonName, 0);
				}

				// Aggregate intensity (like MAG workflow line 2212-2214)
				double[] taxIntensity = taxIntensityMaps[rankIndex].get(taxonName);
				for (int i = 0; i < proteinIntensity.length; i++) {
					taxIntensity[i] += proteinIntensity[i];
				}

				// Aggregate counts
				taxPsmCountMaps[rankIndex].put(taxonName, taxPsmCountMaps[rankIndex].get(taxonName) + psmCount);
				taxPepCountMaps[rankIndex].get(taxonName).addAll(pepCount);
				taxProteinCountMaps[rankIndex].put(taxonName, taxProteinCountMaps[rankIndex].get(taxonName) + 1);
			}
		}

		// Write taxonomy reports (similar structure to MAG workflow)
		try {
			// Write taxa.tsv
			File taxaFile = new File(taxFile, "Taxa.tsv");
			try (PrintWriter writer = new PrintWriter(taxaFile)) {
				StringBuilder header = new StringBuilder();
				header.append("Name\t");
				header.append("Rank");
				for (int i = 0; i < mainRanks.length; i++) {
					header.append("\t").append(mainRanks[i].getName());
				}
				header.append("\tProtein_count\tPSM_count\tPeptide_count");
				for (String intensity : intensityTitle) {
					if (intensity.startsWith("Intensity ")) {
						header.append("\t").append(intensity);
					} else {
						header.append("\tIntensity ").append(intensity);
					}
				}
				writer.println(header);

				// Write data for each rank (like MAG workflow)
				for (int i = 0; i < taxRankMaps.length; i++) {
					for (String taxonName : taxRankMaps[i].keySet()) {
						StringBuilder sb = new StringBuilder();
						sb.append(taxonName).append("\t");
						sb.append(mainRanks[i].getName());

						String[] lineage = taxRankMaps[i].get(taxonName);
						for (int j = 0; j < mainRanks.length; j++) {
							if (j < lineage.length) {
								sb.append("\t").append(lineage[j]);
							} else {
								sb.append("\t");
							}
						}

						sb.append("\t").append(taxProteinCountMaps[i].get(taxonName));
						sb.append("\t").append(taxPsmCountMaps[i].get(taxonName));
						sb.append("\t").append(taxPepCountMaps[i].get(taxonName).size());

						double[] intensity = taxIntensityMaps[i].get(taxonName);
						for (double val : intensity) {
							sb.append("\t").append(val);
						}

						writer.println(sb);
					}
				}
			}

			LOGGER.info(taskName + ": taxonomy analysis completed - wrote " + taxonToProteins.size()
					+ " taxonomic IDs to " + taxaFile.getName());
			System.out.println(format.format(new Date()) + "\t" + taskName + ": taxonomy analysis completed");

		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing taxonomy information", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": error in writing taxonomy information");
			return false;
		}

		return true;
	}

	private boolean generateFinalReports(String taskName) {

		updateMessage(taskName + ": exporting report");

		LOGGER.info(taskName + ": exporting report started");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": exporting report started");

		// Report templates are already copied in the main call() method
		// Just ensure the report directory exists and set the paths
		File reportDir = new File(this.metaPar.getDbSearchResultFile(), "report");
		if (!reportDir.exists()) {
			reportDir.mkdirs();
		}

		File reportHtmlDir = new File(reportDir + "\\reports\\rmd_html");
		if (!reportHtmlDir.exists()) {
			reportHtmlDir.mkdirs();
		}

		// Export metadata
		File summaryMetaFile = new File(this.metaPar.getDbSearchResultFile(), "metainfo.tsv");
		try {
			metaPar.getMetadata().exportMetadata(summaryMetaFile);
		} catch (FileNotFoundException e) {
			LOGGER.error(taskName + ": error in exporting metadata", e);
		}

		// Generate summary report
		MetaPythonReportTask reportTask = new MetaPythonReportTask(threadCount, mrp);
		this.exportSummary(reportTask);

		// Add individual report tasks
		boolean reportSuccess = true;

		// Add summary report
		File summaryReportFile = new File(this.metaPar.getDbSearchResultFile(), "final_summary.tsv");
		if (summaryReportFile.exists() && summaryReportFile.length() > 0) {
			File report_ID_summary = new File(reportHtmlDir, "report_ID_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile,
						report_ID_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.ID_SUMMARY, summaryReportFile,
						report_ID_summary);
			}
		}

		// Add peptide report
		if (final_pep_txt.exists() && final_pep_txt.length() > 0) {
			File report_peptides_summary = new File(reportHtmlDir, "report_peptides_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt,
						report_peptides_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PEPTIDE, final_pep_txt,
						report_peptides_summary);
			}
		}

		// Add protein report
		if (final_pro_txt.exists() && final_pro_txt.length() > 0) {
			File report_proteinGroups_summary = new File(reportHtmlDir, "report_proteinGroups_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
						report_proteinGroups_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.PROTEIN, final_pro_txt,
						report_proteinGroups_summary);
			}
		}

		// Add taxonomy report
		File taxaReportFile = new File(metaPar.getDbSearchResultFile(), "taxonomy_analysis\\Taxa.tsv");
		if (taxaReportFile.exists() && taxaReportFile.length() > 0) {
			File report_taxonomy_summary = new File(reportHtmlDir, "report_taxonomy_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile,
						report_taxonomy_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.TAXONOMY, taxaReportFile,
						report_taxonomy_summary);
			}
		}

		// Add function report
		File funReportTsv = new File(metaPar.getDbSearchResultFile(), "functional_annotation\\functions.tsv");
		if (funReportTsv.exists() && funReportTsv.length() > 0) {
			File report_function_summary = new File(reportHtmlDir, "report_function_summary.html");
			if (summaryMetaFile != null && summaryMetaFile.exists()) {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary, summaryMetaFile);
			} else {
				reportSuccess &= reportTask.addTask(MetaPythonReportTask.FUNCTION, funReportTsv,
						report_function_summary);
			}
		}

		if (!reportSuccess) {
			LOGGER.warn(taskName + ": Some report tasks could not be added");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": Some report tasks could not be added");
		}

		// Run the report generation task
		Thread reportThread = new Thread(reportTask);
		reportThread.start();
		try {
			reportThread.join();
		} catch (InterruptedException e) {
			LOGGER.error(taskName + ": Report generation thread was interrupted", e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": Report generation thread was interrupted");
			Thread.currentThread().interrupt();
			return false;
		}

		try {
			boolean finish = reportTask.get();
			if (finish) {
				LOGGER.info(taskName + ": report generation successes");
				System.out.println(format.format(new Date()) + "\t" + taskName + ": report generation successes");
				return true;
			}
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(taskName + ": error in the report generation", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in the report generation");
			return false;
		}

		return true;
	}
	
	private String buildSpeciesKey(int taxonId) {
		if (taxonId <= 0) {
			return "taxon_unknown";
		}
		return "taxon_" + taxonId;
	}
	


    /**
     * Helper method to update peptide protein mappings with final protein names.
     * 
     * @param peptides Array of peptides to update
     * @param finalProteins Set of final selected protein names
     */
    private void updatePeptideProteinMappings(MetaPeptide[] peptides, HashSet<String> finalProteins) {
        int updatedCount = 0;
        
        for (MetaPeptide peptide : peptides) {
            if (peptide != null && peptide.getProteins() != null) {
                ArrayList<String> newProteins = new ArrayList<String>();
                
                for (String protein : peptide.getProteins()) {
                    if (finalProteins.contains(protein)) {
                        newProteins.add(protein);
                    }
                }
                
                if (!newProteins.isEmpty()) {
                    peptide.setProteins(newProteins.toArray(new String[newProteins.size()]));
                    updatedCount++;
                }
            }
        }
        
        LOGGER.info(getTaskName() + ": Updated " + updatedCount + " peptides with final protein mappings");
    }
}






