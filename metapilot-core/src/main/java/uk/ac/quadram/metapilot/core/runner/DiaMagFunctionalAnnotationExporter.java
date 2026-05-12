package uk.ac.quadram.metapilot.core.runner;

import org.dom4j.DocumentException;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLReader2;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinXMLWriter2;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Headless MAG functional annotation export (TSV + XML), matching the Windows GUI workflow.
 */
public final class DiaMagFunctionalAnnotationExporter {

    private DiaMagFunctionalAnnotationExporter() {
    }

    public static void export(
            File resultDir,
            MagDbItem magDb,
            File funcDbFile,
            String microDbCatalogueId,
            String quantitativeMode,
            MetaProtein[] proteins,
            String[] intensityTitles,
            ProgressCallback cb
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (resultDir == null) throw new IOException("resultDir is null");
        if (magDb == null) throw new IOException("magDb is null");
        if (funcDbFile == null || !funcDbFile.isFile()) {
            throw new IOException("Function definition database not found: "
                    + (funcDbFile == null ? "null" : funcDbFile.getAbsolutePath()));
        }
        if (proteins == null || proteins.length == 0) {
            cb.onMessage("DIA functional annotation: no proteins to annotate");
            return;
        }
        if (microDbCatalogueId == null || microDbCatalogueId.isBlank()) {
            throw new IOException("MAG catalogue id is required for functional annotation");
        }

        File funcDir = new File(resultDir, "functional_annotation");
        if (!funcDir.exists()) funcDir.mkdirs();

        File functionsXml = new File(funcDir, "functions.xml");
        File functionsTsv = new File(funcDir, "functions.tsv");

        if (!functionsXml.exists() || functionsXml.length() == 0) {
            cb.onMessage("DIA functional annotation: annotating proteins");
            annotateProteins(magDb, funcDbFile, microDbCatalogueId, quantitativeMode, proteins, intensityTitles, functionsXml);
        } else {
            cb.onMessage("DIA functional annotation: reusing existing " + functionsXml.getName());
        }

        if (!functionsXml.exists() || functionsXml.length() == 0) {
            throw new IOException("Functional annotation did not produce: " + functionsXml.getAbsolutePath());
        }

        MetaProteinXMLReader2 reader = new MetaProteinXMLReader2(functionsXml);
        if (!functionsTsv.exists()) {
            cb.onMessage("DIA functional annotation: writing " + functionsTsv.getName());
            reader.exportTsv(functionsTsv);
            reader.exportGOReport(new File(funcDir, "Gene_Ontology_id_report.tsv"));
            reader.exportECReport(new File(funcDir, "EC_id_report.tsv"));
            reader.exportKeggReport(new File(funcDir, "KEGG_ko_report.tsv"));
            reader.exportBriteReport(new File(funcDir, "BRITE_report.tsv"));
            reader.exportCOGReport(new File(funcDir, "COG_accession_report.tsv"));
            reader.exportNOGReport(new File(funcDir, "NOG_accession_report.tsv"));
            reader.exportCOGCategoryReport(new File(funcDir, "COG_category_report.tsv"));
            reader.exportNOGCategoryReport(new File(funcDir, "NOG_category_report.tsv"));
        } else {
            cb.onMessage("DIA functional annotation: reusing existing " + functionsTsv.getName());
        }
    }

    private static void annotateProteins(
            MagDbItem magDb,
            File funcDbFile,
            String microDbCatalogueId,
            String quantitativeMode,
            MetaProtein[] proteins,
            String[] intensityTitles,
            File functionsXml
    ) throws SQLException, IOException, DocumentException {
        HashMap<String, ArrayList<MetaProtein>> genomeProMap = new HashMap<>();
        for (MetaProtein pro : proteins) {
            if (pro == null) continue;
            String name = pro.getName();
            if (name == null || !name.startsWith(magDb.getIdentifier())) continue;
            String genome = name.split("_")[0];
            if (genome.endsWith(".1")) {
                genome = genome.substring(0, genome.length() - 2);
            }
            genomeProMap.computeIfAbsent(genome, k -> new ArrayList<>()).add(pro);
        }

        HashMap<String, String[]> usedCogMap = new HashMap<>();
        HashMap<String, String[]> usedNogMap = new HashMap<>();
        HashMap<String, String> usedKeggMap = new HashMap<>();
        HashMap<String, uk.ac.quadram.metapilot.core.function.v2.GoObo> usedGoMap = new HashMap<>();
        HashMap<String, uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission> usedEcMap = new HashMap<>();
        ArrayList<MetaProteinAnnoMag> annotated = new ArrayList<>();

        for (String genome : genomeProMap.keySet()) {
            ArrayList<MetaProtein> genomeProteins = genomeProMap.get(genome);
            MetaProtein[] batch = genomeProteins.toArray(new MetaProtein[0]);
            MagFuncSearcher searcher = new MagFuncSearcher(
                    microDbCatalogueId, funcDbFile.getAbsolutePath(), magDb, genome);
            MetaProteinAnnoMag[] matches = searcher.match(batch);
            if (matches != null) {
                for (MetaProteinAnnoMag match : matches) {
                    annotated.add(match);
                }
                usedCogMap.putAll(searcher.getUsedCogMap());
                usedNogMap.putAll(searcher.getUsedNogMap());
                usedKeggMap.putAll(searcher.getUsedKeggMap());
                usedEcMap.putAll(searcher.getUsedEcMap());
                usedGoMap.putAll(searcher.getUsedGoMap());
            }
        }

        String quantMode = quantitativeMode == null || quantitativeMode.isBlank() ? "Label-free" : quantitativeMode;
        MetaProteinXMLWriter2 writer = new MetaProteinXMLWriter2(
                functionsXml.getAbsolutePath(), MetaConstants.pFindMag, quantMode, intensityTitles);
        writer.addProteins(
                annotated.toArray(new MetaProteinAnnoMag[0]),
                usedCogMap,
                usedNogMap,
                usedKeggMap,
                usedGoMap,
                usedEcMap);
        writer.close();
    }
}
