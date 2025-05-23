package uk.ac.quadram.metapilot.task.mag.gui;

import java.awt.EventQueue;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;

import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.task.hgm.gui.MetaPilotMainFrameHGM;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.par.MetaParaIOMag;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.mag.par.MetaSourcesIoMag;
import uk.ac.quadram.metapilot.task.mag.par.MetaSourcesMag;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.pfind.gui.MetaPilotLicenseFrame;
import uk.ac.quadram.metapilot.task.v2.gui.MetaLabSourcePathDialog;
import net.miginfocom.swing.MigLayout;
import uk.ac.quadram.metapilot.task.MetaPilotTask;

public class MetaPilotMainFrameMag extends MetaPilotMainFrameHGM {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3053879427444409502L;
	protected static Logger LOGGER = LogManager.getLogger(MetaPilotMainFrameMag.class);
	private MetaPilotMainPanelMag contentPane;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {

		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException | InstantiationException | IllegalAccessException
				| UnsupportedLookAndFeelException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		EventQueue.invokeLater(new Runnable() {
			public void run() {
				MetaPilotMainFrameMag frame = null;
				try {
					frame = new MetaPilotMainFrameMag();
					frame.setVisible(true);
				} catch (Exception e) {
					MetaPilotMainFrameMag.LOGGER.error(e);
				}
			}
		});
	}

	/**
	 * Create the frame.
	 * 
	 * @throws URISyntaxException
	 * @throws IOException
	 */
	public MetaPilotMainFrameMag() throws URISyntaxException, IOException {
		super();
	}

	protected void initial() throws IOException, URISyntaxException {

		File parameter = new File(MetaPilotTask.getApplicationBaseDirectory(), "parameter" + MetaParaIOMag.versionFile);

		File resources = new File(MetaPilotTask.getApplicationBaseDirectory(), "resources" + MetaParaIOMag.versionFile);

		File dbpar = new File(MetaPilotTask.getApplicationBaseDirectory(), "dbpar" + MetaParaIOMag.versionFile);
		
		MetaSourcesMag msv2;
		if (resources.exists()) {
			msv2 = MetaSourcesIoMag.parse(resources);
		} else {
			File oldResources = new File(MetaPilotTask.getApplicationBaseDirectory(), "resources" + MetaParaIOMag.oldVersionFiles[0]);

			if (oldResources.exists()) {
				msv2 = MetaSourcesIoMag.parse(oldResources);
				msv2.setVersion(MetaParaIOMag.versionFile);
				MetaSourcesIoMag.export(msv2, resources);
			} else {
				msv2 = MetaSourcesIoMag.parse(resources);
			}
		}

		MetaBaseParameter parV3 = MetaParaIOMag.parse(parameter);

		MagDbItem[][] magDbItems = MagDbConfigIO.parse(dbpar);

		((MetaParameterMag) parV3).setMagDbItems(magDbItems);

		getContentPane().setLayout(new MigLayout("", "[1200:1600:1920,grow 600]", "[800:900:1080,grow]"));

		LOGGER.info("Initial the content panel...");
		
		this.contentPane = new MetaPilotMainPanelMag(parV3, (new MetaResourceIO()).loadParameters());

		getContentPane().add(contentPane, "cell 0 0,grow");

		setIconImage(Toolkit.getDefaultToolkit()
				.getImage(MetaPilotMainFrameMag.class.getResource("/MetaPilot.png")));

		setTitle("MetaLab " + MetaParaIOMag.version);

		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);

		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);

		JMenuItem mntmLoadParameter = new JMenuItem("Load parameter");
		mntmLoadParameter.setIcon(
				new ImageIcon(MetaPilotMainFrameMag.class.getResource("/toolbarButtonGraphics/general/Open16.gif")));

		mntmLoadParameter.addActionListener(l -> {
			JFileChooser fileChooser = new JFileChooser();
			fileChooser.setDialogTitle("Load parameter file");
			fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			fileChooser.addChoosableFileFilter(new FileFilter() {
				public String getDescription() {
					return "MetaLab parameter file (.json)";
				}

				public boolean accept(File f) {
					if (f.getName().endsWith("json")) {
						return true;
					} else {
						return false;
					}
				}
			});

			int returnValue = fileChooser.showOpenDialog(MetaPilotMainFrameMag.this);
			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = fileChooser.getSelectedFile();
				MetaParameterMag loadPar = MetaParaIOMag.parse(file);
				if (loadPar.getWorkflowType() == MetaLabWorkflowType.pFindMAG) {

					MetaPilotMainFrameMag.this.getContentPane().remove(this.contentPane);
					contentPane = new MetaPilotMainPanelMag(loadPar, (new MetaResourceIO()).loadParameters());
					MetaPilotMainFrameMag.this.getContentPane().add(contentPane, "cell 0 0,grow");
					MetaPilotMainFrameMag.this.validate();

				} else {
					JOptionPane.showMessageDialog(this,
							"Only parameters with workflow type " + MetaLabWorkflowType.pFindMAG.getDescription()
									+ " is supported, the workflow type you loaded is "
									+ loadPar.getWorkflowType().getDescription(),
							"Warning", JOptionPane.WARNING_MESSAGE);
				}
			}
		});
		mnFile.add(mntmLoadParameter);

		JMenuItem mntmSaveParameter = new JMenuItem("Save parameter");
		mntmSaveParameter.setIcon(
				new ImageIcon(MetaPilotMainFrameMag.class.getResource("/toolbarButtonGraphics/general/Save16.gif")));
		mnFile.add(mntmSaveParameter);
		mntmSaveParameter.addActionListener(l -> {
			JFileChooser fileChooser = new JFileChooser();
			fileChooser.setDialogTitle("Save parameter file");

			int value = fileChooser.showSaveDialog(MetaPilotMainFrameMag.this);

			if (value == JFileChooser.APPROVE_OPTION) {
				File fileToSave = fileChooser.getSelectedFile();
				if (!fileToSave.getName().endsWith(".json")) {
					fileToSave = new File(fileToSave.getParent(), fileToSave.getName() + ".json");
				}

				((MetaPilotMainPanelMag) contentPane).updateParameter();
				MetaBaseParameter par = ((MetaPilotMainPanelMag) contentPane).getParameter();

				MetaParaIOMag.export((MetaParameterMag) par, fileToSave);
			}
		});

		JMenuItem mntmExit = new JMenuItem("Exit");
		mntmExit.setIcon(
				new ImageIcon(MetaPilotMainFrameMag.class.getResource("/toolbarButtonGraphics/general/Stop16.gif")));
		mntmExit.addActionListener(l -> {
			this.dispose();
		});
		mnFile.add(mntmExit);

		JMenu mnTools = new JMenu("Tools");
		menuBar.add(mnTools);

		JMenuItem mntmVerify = new JMenuItem("Activation");
		mntmVerify.setIcon(
				new ImageIcon(MetaPilotMainFrameMag.class.getResource("/toolbarButtonGraphics/general/Add16.gif")));
		mntmVerify.addActionListener(l -> {
			MetaPilotLicenseFrame frame;
			try {
				frame = new MetaPilotLicenseFrame();
				frame.setAlwaysOnTop(true);
				frame.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
				frame.setVisible(true);
			} catch (URISyntaxException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		});
		mnTools.add(mntmVerify);

		JMenu mnSetting = new JMenu("Setting");
		menuBar.add(mnSetting);

		JMenuItem mntmResource = new JMenuItem("Resource");
		mntmResource.setIcon(
				new ImageIcon(MetaPilotMainFrameMag.class.getResource("/toolbarButtonGraphics/general/Bookmarks16.gif")));
		mntmResource.addActionListener(l -> {
			MetaLabSourcePathDialog dialog = new MetaLabSourceDialogMag(msv2, resources, this);
			dialog.setAlwaysOnTop(true);
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		});
		mnSetting.add(mntmResource);

		JMenu mnHelp = new JMenu("Help");
		menuBar.add(mnHelp);

		JMenuItem mntmHelp = new JMenuItem("Help");
		mntmHelp.setIcon(
				new ImageIcon(MetaPilotMainFrameMag.class.getResource("/toolbarButtonGraphics/general/Help16.gif")));
		mntmHelp.addActionListener(l -> {
			try {
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("https://wiki.imetalab.ca/"));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		});
		mnHelp.add(mntmHelp);

		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {

				((MetaPilotMainPanelMag) contentPane).updateParameter();
				MetaBaseParameter par = ((MetaPilotMainPanelMag) contentPane).getParameter();

				MetaParaIOMag.export((MetaParameterMag) par, parameter);

				MagDbItem[][] magDbItems = ((MetaParameterMag) par).getMagDbItems();

				MagDbConfigIO.export(magDbItems, dbpar);

				e.getWindow().dispose();
			}
		});
		
		LOGGER.info("Finished...");
	}

}
