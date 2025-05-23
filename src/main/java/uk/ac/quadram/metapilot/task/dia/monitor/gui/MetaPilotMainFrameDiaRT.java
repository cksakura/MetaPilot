package uk.ac.quadram.metapilot.task.dia.monitor.gui;

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
import javax.swing.JPanel;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;

import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.task.dia.gui.MetaPilotMainFrameDia;
import uk.ac.quadram.metapilot.task.dia.gui.MetaLabSourceDialogDia;
import uk.ac.quadram.metapilot.task.dia.monitor.par.MetaParIODiaRT;
import uk.ac.quadram.metapilot.task.dia.monitor.par.MetaParameterDiaRT;
import uk.ac.quadram.metapilot.task.dia.par.MetaParaIoDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaSourcesDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaSourcesIoDia;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.pfind.gui.MetaPilotLicenseFrame;
import uk.ac.quadram.metapilot.task.v2.gui.MetaLabSourcePathDialog;
import net.miginfocom.swing.MigLayout;

public class MetaPilotMainFrameDiaRT extends MetaPilotMainFrameDia {

	/**
	 * 
	 */
	private static final long serialVersionUID = 8710147856572038547L;
	private JPanel contentPane;

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
				try {
					MetaPilotMainFrameDiaRT frame = new MetaPilotMainFrameDiaRT();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
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
	public MetaPilotMainFrameDiaRT() throws URISyntaxException, IOException {
		super();
	}

	protected void initial() throws IOException, URISyntaxException {

		File parameter = new File(ClassLoader.getSystemClassLoader().getResource(".").toURI().getPath(),
				"parameters" + MetaParIODiaRT.versionFile);

		File resources = new File(ClassLoader.getSystemClassLoader().getResource(".").toURI().getPath(),
				"resources" + MetaParaIoDia.versionFile);

		File dbpar = new File(ClassLoader.getSystemClassLoader().getResource(".").toURI().getPath(),
				MagDbConfigIO.versionFile);

		MetaSourcesDia msv = MetaSourcesIoDia.parse(resources);
		MagDbItem[][] magDbItems = MagDbConfigIO.parse(dbpar);
		MetaParameterDiaRT par = MetaParIODiaRT.parse(parameter);
		par.setMagDbItems(magDbItems);

		getContentPane().setLayout(new MigLayout("", "[1200:1600:1920,grow]", "[800:900:1080,grow]"));

		this.contentPane = new MetaPilotMainPanelDiaRT(par, (new MetaResourceIO()).loadParameters());

		getContentPane().add(contentPane, "cell 0 0,grow");

		setIconImage(Toolkit.getDefaultToolkit()
				.getImage(MetaPilotMainFrameDiaRT.class.getResource("/MetaPilot.png")));

		setTitle("MetaLab " + MetaParIODiaRT.version);

		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);

		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);

		JMenuItem mntmLoadParameter = new JMenuItem("Load parameter");
		mntmLoadParameter.setIcon(
				new ImageIcon(MetaPilotMainFrameDia.class.getResource("/toolbarButtonGraphics/general/Open16.gif")));

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

			int returnValue = fileChooser.showOpenDialog(MetaPilotMainFrameDiaRT.this);
			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = fileChooser.getSelectedFile();
				MetaParameterDiaRT loadPar = MetaParIODiaRT.parse(file);
				if (loadPar.getWorkflowType() == MetaLabWorkflowType.DiaNNMAGRT) {
					MetaPilotMainFrameDiaRT.this.getContentPane().remove(this.contentPane);
					contentPane = new MetaPilotMainPanelDiaRT(loadPar, (new MetaResourceIO()).loadParameters());
					MetaPilotMainFrameDiaRT.this.getContentPane().add(contentPane, "cell 0 0,grow");
					MetaPilotMainFrameDiaRT.this.validate();

				} else {
					JOptionPane.showMessageDialog(this,
							"Only parameters with workflow type " + MetaLabWorkflowType.DiaNNMAGRT.getDescription()
									+ " is supported, the workflow type you loaded is "
									+ loadPar.getWorkflowType().getDescription(),
							"Warning", JOptionPane.WARNING_MESSAGE);
				}
			}
		});
		mnFile.add(mntmLoadParameter);

		JMenuItem mntmSaveParameter = new JMenuItem("Save parameter");
		mntmSaveParameter.setIcon(
				new ImageIcon(MetaPilotMainFrameDia.class.getResource("/toolbarButtonGraphics/general/Save16.gif")));
		mnFile.add(mntmSaveParameter);
		mntmSaveParameter.addActionListener(l -> {
			JFileChooser fileChooser = new JFileChooser();
			fileChooser.setDialogTitle("Save parameter file");

			int value = fileChooser.showSaveDialog(MetaPilotMainFrameDiaRT.this);

			if (value == JFileChooser.APPROVE_OPTION) {
				File fileToSave = fileChooser.getSelectedFile();
				if (!fileToSave.getName().endsWith(".json")) {
					fileToSave = new File(fileToSave.getParent(), fileToSave.getName() + ".json");
				}

				((MetaPilotMainPanelDiaRT) contentPane).updateParameter();
				MetaBaseParameter savePar = ((MetaPilotMainPanelDiaRT) contentPane).getParameter();

				MetaParIODiaRT.export((MetaParameterDiaRT) savePar, fileToSave);
			}
		});

		JMenuItem mntmExit = new JMenuItem("Exit");
		mntmExit.setIcon(
				new ImageIcon(MetaPilotMainFrameDia.class.getResource("/toolbarButtonGraphics/general/Stop16.gif")));
		mntmExit.addActionListener(l -> {
			this.dispose();
		});
		mnFile.add(mntmExit);

		JMenu mnTools = new JMenu("Tools");
		menuBar.add(mnTools);

		JMenuItem mntmVerify = new JMenuItem("Activation");
		mntmVerify.setIcon(
				new ImageIcon(MetaPilotMainFrameDia.class.getResource("/toolbarButtonGraphics/general/Add16.gif")));
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
				new ImageIcon(MetaPilotMainFrameDia.class.getResource("/toolbarButtonGraphics/general/Bookmarks16.gif")));
		mntmResource.addActionListener(l -> {
			MetaLabSourcePathDialog dialog = new MetaLabSourceDialogDia((MetaSourcesDia) msv, resources, this);
			dialog.setAlwaysOnTop(true);
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		});
		mnSetting.add(mntmResource);

		JMenu mnHelp = new JMenu("Help");
		menuBar.add(mnHelp);

		JMenuItem mntmHelp = new JMenuItem("Help");
		mntmHelp.setIcon(
				new ImageIcon(MetaPilotMainFrameDiaRT.class.getResource("/toolbarButtonGraphics/general/Help16.gif")));
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

				((MetaPilotMainPanelDiaRT) contentPane).updateParameter();
				MetaBaseParameter par = ((MetaPilotMainPanelDiaRT) contentPane).getParameter();

				MetaParIODiaRT.export((MetaParameterDiaRT) par, parameter);

				MagDbItem[][] magDbItems = ((MetaParameterDiaRT) par).getMagDbItems();

				MagDbConfigIO.export(magDbItems, dbpar);

				e.getWindow().dispose();
			}
		});
	}

}
