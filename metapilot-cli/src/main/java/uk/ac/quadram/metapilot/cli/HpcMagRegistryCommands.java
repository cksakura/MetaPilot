package uk.ac.quadram.metapilot.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.File;
import java.nio.file.Path;
import java.util.Date;

final class HpcMagRegistryCommands {

    private static final String MGNIFY_DB_FOLDER = "original_db";
    private static final String MGNIFY_FUNC_FOLDER = "eggNOG";
    private static final String MGNIFY_TAXA_FILE = "genomes-all_metadata.tsv";

    private HpcMagRegistryCommands() {
    }

    static int runMagList(String[] args) {
        Options options = buildMagListOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printMagListHelp(options);
            return 2;
        }
        if (cmd.hasOption("help")) {
            printMagListHelp(options);
            return 0;
        }

        try {
            HpcSiteConfig.ensureSiteConfig(false);
            Path magdbCfg = HpcSiteConfig.resolveMagDbConfigPath(cmd.getOptionValue("magdb"));
            MagDbItem[][] parsed = MagDbConfigIO.parse(magdbCfg.toFile());
            if (parsed == null) {
                System.err.println("ERROR: failed to read MAG registry: " + magdbCfg);
                return 1;
            }
            System.out.println("MAG registry: " + magdbCfg);
            printGroup("MGnify", parsed[0]);
            printGroup("localMAGs", parsed[1]);
            return 0;
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return 1;
        }
    }

    static int runMagRegister(String[] args) {
        Options options = buildMagRegisterOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printMagRegisterHelp(options);
            return 2;
        }
        if (cmd.hasOption("help")) {
            printMagRegisterHelp(options);
            return 0;
        }

        String catalogue = cmd.getOptionValue("catalogue");
        String version = cmd.getOptionValue("version");
        String repository = cmd.getOptionValue("repository");
        if (catalogue == null || version == null || repository == null) {
            System.err.println("ERROR: --catalogue, --version, and --repository are required.");
            printMagRegisterHelp(options);
            return 2;
        }

        String identifier = cmd.getOptionValue("identifier", "MGYG");
        boolean mgnify = !cmd.hasOption("local");

        try {
            HpcSiteConfig.ensureSiteConfig(false);
            Path magdbCfg = HpcSiteConfig.resolveMagDbConfigPath(cmd.getOptionValue("magdb"));
            File repo = new File(repository.trim());
            if (!repo.isDirectory()) {
                System.err.println("ERROR: repository is not a directory: " + repo.getAbsolutePath());
                return 2;
            }

            MagDbItem[][] parsed = MagDbConfigIO.parse(magdbCfg.toFile());
            if (parsed == null) {
                System.err.println("ERROR: failed to read MAG registry: " + magdbCfg);
                return 1;
            }

            MagDbItem item = findItem(parsed, mgnify, catalogue.trim(), identifier.trim());
            boolean saved;
            if (item == null) {
                item = new MagDbItem(catalogue.trim(), identifier.trim(), repo);
                addRegisteredVersion(item, version.trim());
                saved = MagDbConfigIO.addMagDbItem(item, mgnify, magdbCfg.toFile());
            } else {
                item.setRepositoryFile(repo);
                boolean versionFound = false;
                for (MagDbItem.MagVersion magVersion : item.getVersionList()) {
                    if (version.trim().equals(magVersion.getCatalogueVersion())) {
                        item.setUsedVersion(magVersion);
                        item.setAvailable();
                        versionFound = true;
                        break;
                    }
                }
                if (!versionFound) {
                    addRegisteredVersion(item, version.trim());
                }
                saved = MagDbConfigIO.updateMagDbItem(item, mgnify, magdbCfg.toFile());
            }

            if (!saved) {
                System.err.println("ERROR: failed to update MAG registry: " + magdbCfg);
                return 1;
            }

            System.out.println("Registered MAG catalogue " + catalogue + " / " + version + " at " + repo.getAbsolutePath());
            return 0;
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static void addRegisteredVersion(MagDbItem item, String version) {
        item.addVersion(version, 0, new Date(), true, MGNIFY_DB_FOLDER, MGNIFY_FUNC_FOLDER, MGNIFY_TAXA_FILE);
        for (MagDbItem.MagVersion magVersion : item.getVersionList()) {
            if (version.equals(magVersion.getCatalogueVersion())) {
                item.setUsedVersion(magVersion);
                break;
            }
        }
    }

    private static MagDbItem findItem(MagDbItem[][] parsed, boolean mgnify, String catalogue, String identifier) {
        MagDbItem[] group = mgnify ? parsed[0] : parsed[1];
        if (group == null) {
            return null;
        }
        for (MagDbItem item : group) {
            if (item != null && catalogue.equals(item.getCatalogueID()) && identifier.equals(item.getIdentifier())) {
                return item;
            }
        }
        return null;
    }

    private static void printGroup(String label, MagDbItem[] items) {
        System.out.println(label + ":");
        if (items == null || items.length == 0) {
            System.out.println("  (none)");
            return;
        }
        for (MagDbItem item : items) {
            if (item == null) {
                continue;
            }
            File repo = item.getRepositoryFile();
            System.out.println("  " + item.getCatalogueID() + " [" + item.getIdentifier() + "]"
                    + (repo != null ? " repo=" + repo.getAbsolutePath() : ""));
            for (MagDbItem.MagVersion version : item.getVersionList()) {
                System.out.println("    - " + version.getCatalogueVersion()
                        + (version.isLocalAvailable() ? " (local)" : " (remote)"));
            }
        }
    }

    private static Options buildMagListOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h").longOpt("help").desc("Show help").build());
        opts.addOption(Option.builder().longOpt("magdb").hasArg().argName("FILE").desc("MAG registry JSON override").build());
        return opts;
    }

    private static Options buildMagRegisterOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h").longOpt("help").desc("Show help").build());
        opts.addOption(Option.builder().longOpt("catalogue").hasArg().argName("ID").desc("Catalogue ID").required().build());
        opts.addOption(Option.builder().longOpt("version").hasArg().argName("VER").desc("Catalogue version").required().build());
        opts.addOption(Option.builder().longOpt("repository").hasArg().argName("DIR").desc("Parent directory containing catalogue/version folders").required().build());
        opts.addOption(Option.builder().longOpt("identifier").hasArg().argName("ID").desc("Catalogue identifier (default: MGYG)").build());
        opts.addOption(Option.builder().longOpt("local").desc("Register under localMAGs instead of MGnify").build());
        opts.addOption(Option.builder().longOpt("magdb").hasArg().argName("FILE").desc("MAG registry JSON override").build());
        return opts;
    }

    private static void printMagListHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp("java -jar metapilot-cli-<version>-shaded.jar mag-list [options]",
                "\nList catalogues and locally available versions from the site MAG registry.",
                options, "", true);
    }

    private static void printMagRegisterHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp("java -jar metapilot-cli-<version>-shaded.jar mag-register [options]",
                "\nRegister or update a local MAG catalogue repository path in the site MAG registry.",
                options, "", true);
    }
}
