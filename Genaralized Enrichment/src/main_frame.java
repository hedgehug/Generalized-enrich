import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JButton;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.awt.event.ActionEvent;
import java.awt.Font;
import javax.swing.JTextField;
import javax.swing.JTabbedPane;
import javax.swing.JLabel;
import javax.swing.JFileChooser;
import javax.swing.JCheckBox;
import javax.swing.SwingConstants;


public class main_frame 
{

	private JFrame frame;
	private JTextField experimentFilePath;
	private JTextField targetFilePath;
	private JTextField experiment_identifier_1;
	private JTextField experiment_identifier_2;
	private JTextField target_identifier_1;
	private JTextField target_identifier_2;
	
	private String experiment_file_path;
	private String target_file_path;
	static String work_directory;
	private ExpressionFileLoader experiment_expression_file_loader;
	private ExpressionFileLoader target_expression_file_loader;
	static int num_permutation;
	public static boolean is_gene_list;
	private GeneExpression experiment_gene_expression;
	private GeneExpression target_gene_expression;
	private GeneListFileLoader target_gene_list;
	static int window_stride;
	static int initial_window_start;
	private Comparator comparator;
	
	private String experiment_sample_identifier_1;
	private String experiment_sample_identifier_2;
	private String target_sample_identifier_1;
	private String target_sample_identifier_2;
	
	private double max_ratio;
	private double significance;
	
	private JTextField permutation_number;
	private JTextField working_directory;
	private JTextField txtWindowSize;
	private JTextField txtInitialWindowSize;
	
	

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					main_frame window = new main_frame();
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public main_frame() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		is_gene_list = false;
		num_permutation = 100;
		window_stride = 10;
		initial_window_start = 100;
		
		frame = new JFrame();
		frame.setBounds(100, 100, 921, 630);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setLayout(null);
		
		JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP);
		tabbedPane.setBounds(0, 0, 921, 608);
		frame.getContentPane().add(tabbedPane);
		
		JPanel tabLoadData = new JPanel();
		tabbedPane.addTab("Load Data", null, tabLoadData, null);
		tabLoadData.setLayout(null);
		
		JLabel lblLoadExperiment = new JLabel("Load Experiment");
		lblLoadExperiment.setFont(new Font("Lucida Grande", Font.PLAIN, 20));
		lblLoadExperiment.setBounds(6, 6, 169, 32);
		tabLoadData.add(lblLoadExperiment);
		
		experimentFilePath = new JTextField();
		experimentFilePath.setBounds(187, 6, 543, 32);
		tabLoadData.add(experimentFilePath);
		experimentFilePath.setColumns(10);
		
		
		
		JButton btnChooseExperimentFile = new JButton("Browse");
		btnChooseExperimentFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				// open file chooser
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.showOpenDialog(null);
				File experimentFile = fileChooser.getSelectedFile();
				experimentFilePath.setText(experimentFile.getAbsolutePath());		
			}
		});
		
		btnChooseExperimentFile.setBounds(742, 6, 117, 32);
		tabLoadData.add(btnChooseExperimentFile);
		
		JLabel lblNewLabel = new JLabel("Load Target");
		lblNewLabel.setFont(new Font("Lucida Grande", Font.PLAIN, 20));
		lblNewLabel.setBounds(6, 123, 169, 32);
		tabLoadData.add(lblNewLabel);
		
		targetFilePath = new JTextField();
		targetFilePath.setBounds(187, 126, 543, 33);
		tabLoadData.add(targetFilePath);
		targetFilePath.setColumns(10);
		
		JButton btnChooseTargetFile = new JButton("Browse");
		btnChooseTargetFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				// open file chooser
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.showOpenDialog(null);
				File targetFile = fileChooser.getSelectedFile();
				targetFilePath.setText(targetFile.getAbsolutePath());		
			}
		});
		btnChooseTargetFile.setBounds(742, 129, 117, 29);
		tabLoadData.add(btnChooseTargetFile);
		
		JLabel lblId_2 = new JLabel("id 1:");
		lblId_2.setBounds(197, 50, 34, 16);
		tabLoadData.add(lblId_2);
		
		JLabel lblId_1 = new JLabel("id 2:");
		lblId_1.setBounds(197, 76, 34, 16);
		tabLoadData.add(lblId_1);
		
		experiment_identifier_1 = new JTextField();
		experiment_identifier_1.setToolTipText("e.g. prefix for control");
		experiment_identifier_1.setBounds(243, 45, 130, 26);
		tabLoadData.add(experiment_identifier_1);
		experiment_identifier_1.setColumns(10);
		
		experiment_identifier_2 = new JTextField();
		experiment_identifier_2.setToolTipText("e.g. prefix for treatment");
		experiment_identifier_2.setColumns(10);
		experiment_identifier_2.setBounds(243, 71, 130, 26);
		tabLoadData.add(experiment_identifier_2);
		
		JLabel label = new JLabel("id 1:");
		label.setBounds(197, 171, 34, 16);
		tabLoadData.add(label);
		
		JLabel label_1 = new JLabel("id 2:");
		label_1.setBounds(197, 199, 34, 16);
		tabLoadData.add(label_1);
		
		target_identifier_1 = new JTextField();
		target_identifier_1.setToolTipText("e.g. prefix for control");
		target_identifier_1.setColumns(10);
		target_identifier_1.setBounds(243, 166, 130, 26);
		tabLoadData.add(target_identifier_1);
		
		target_identifier_2 = new JTextField();
		target_identifier_2.setToolTipText("e.g. prefix for treatment");
		target_identifier_2.setColumns(10);
		target_identifier_2.setBounds(243, 194, 130, 26);
		tabLoadData.add(target_identifier_2);
		
		JCheckBox chckbxGeneList = new JCheckBox("Gene List");
		chckbxGeneList.setBounds(187, 228, 128, 23);
		tabLoadData.add(chckbxGeneList);
		
		JPanel tabRunAnalysis = new JPanel();
		tabbedPane.addTab("Run Analysis", null, tabRunAnalysis, null);
		tabRunAnalysis.setLayout(null);
		
		JLabel lblRequiredField = new JLabel("Required Field");
		lblRequiredField.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblRequiredField.setBounds(6, 6, 135, 20);
		tabRunAnalysis.add(lblRequiredField);
		
		JLabel lblNewLabel_1 = new JLabel("Permutation");
		lblNewLabel_1.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblNewLabel_1.setBounds(6, 50, 135, 22);
		tabRunAnalysis.add(lblNewLabel_1);
		
		permutation_number = new JTextField();
		permutation_number.setText("1000");
		permutation_number.setBounds(165, 50, 201, 22);
		tabRunAnalysis.add(permutation_number);
		permutation_number.setColumns(10);
		
		JPanel tabSeeResults = new JPanel();
		tabbedPane.addTab("See Results", null, tabSeeResults, null);
		tabSeeResults.setLayout(null);
		
		JLabel lblMaxOddsRatio = new JLabel("Max odds ratio:");
		lblMaxOddsRatio.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblMaxOddsRatio.setBounds(6, 6, 529, 20);
		tabSeeResults.add(lblMaxOddsRatio);
		
		JLabel lblEnrichmentSignificance = new JLabel("Enrichment significance:");
		lblEnrichmentSignificance.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblEnrichmentSignificance.setBounds(6, 36, 529, 20);
		tabSeeResults.add(lblEnrichmentSignificance);
		
		JLabel lblResultsSavedIn = new JLabel("Results saved in:");
		lblResultsSavedIn.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblResultsSavedIn.setBounds(6, 66, 529, 20);
		tabSeeResults.add(lblResultsSavedIn);
		
		JButton btnOpen = new JButton("Open");
		btnOpen.addActionListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent e) 
			{
				// TODO remain to be fixed
				String temp_work_dir = "";
				if(work_directory.contains(" "))
				{
					temp_work_dir = work_directory.replace(" ", "\\ ");
				}
				String command = "/usr/bin/open "+temp_work_dir;
				ProcessBuilder processBuilder = new ProcessBuilder();
				processBuilder.command(command);
				try {
					Process process = processBuilder.start();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				// process.waitFor();
			}
		});
		btnOpen.setBounds(543, 64, 117, 29);
		tabSeeResults.add(btnOpen);
		
		JLabel lblLoadingData = new JLabel("");
		lblLoadingData.setHorizontalAlignment(SwingConstants.CENTER);
		lblLoadingData.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblLoadingData.setBounds(216, 422, 470, 32);
		tabLoadData.add(lblLoadingData);	
		
		
		JButton btnLoadData = new JButton("Load Data");
		btnLoadData.addActionListener(new ActionListener() {
			
			public void actionPerformed(ActionEvent e) {
				
				experiment_file_path = experimentFilePath.getText();
				target_file_path = targetFilePath.getText();
							
				experiment_sample_identifier_1 = experiment_identifier_1.getText();
				experiment_sample_identifier_2 = experiment_identifier_2.getText();
				
				if (chckbxGeneList.isSelected())
				{
					is_gene_list = true;
				}
				else
				{
					is_gene_list = false;
				}
				
				try {
					lblLoadingData.setText("Loading data...");	
					boolean files_loaded = ExpressionLoader();
					if (!files_loaded)
					{
						lblLoadingData.setText("File cannot be loaded! Please check file format.");
					}
					else
					{
						lblLoadingData.setText("Files loaded!");
					}
				} catch (IOException e2) {
					// TODO Auto-generated catch block
					e2.printStackTrace();
				}	
			}
		});
		btnLoadData.setBounds(387, 380, 117, 29);
		tabLoadData.add(btnLoadData);
		
		
		
		JLabel lblDirectory = new JLabel("Directory");
		lblDirectory.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblDirectory.setBounds(6, 80, 84, 20);
		tabRunAnalysis.add(lblDirectory);
		
		working_directory = new JTextField();
		working_directory.setBounds(165, 80, 201, 22);
		tabRunAnalysis.add(working_directory);
		working_directory.setColumns(10);
		
		JLabel lblWindowSize = new JLabel("Window Stride");
		lblWindowSize.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblWindowSize.setBounds(6, 110, 140, 20);
		tabRunAnalysis.add(lblWindowSize);
		
		txtWindowSize = new JTextField();
		txtWindowSize.setText("10");
		txtWindowSize.setBounds(165, 110, 201, 22);
		tabRunAnalysis.add(txtWindowSize);
		txtWindowSize.setColumns(10);
		
		JLabel lblInitialWindowSize = new JLabel("Initial Window Size");
		lblInitialWindowSize.setFont(new Font("Lucida Grande", Font.PLAIN, 16));
		lblInitialWindowSize.setBounds(6, 140, 153, 20);
		tabRunAnalysis.add(lblInitialWindowSize);
		
		txtInitialWindowSize = new JTextField();
		txtInitialWindowSize.setText("100");
		txtInitialWindowSize.setBounds(165, 140, 201, 22);
		tabRunAnalysis.add(txtInitialWindowSize);
		txtInitialWindowSize.setColumns(10);
		
		
		JButton Run_Analysis = new JButton("Run Analysis");
		Run_Analysis.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				work_directory = working_directory.getText();
				num_permutation = Integer.valueOf(permutation_number.getText());
				File dir = new File(work_directory);
				dir.mkdir();
				
				window_stride = Integer.valueOf(txtWindowSize.getText());
				initial_window_start = Integer.valueOf(txtInitialWindowSize.getText());
				
				// calculate statistics for expression				
				CalculateExpression();
				
				// compare gene expression
				try {
					Compare();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				
				lblMaxOddsRatio.setText("Max odds ratio: "+String.format("%.3f", max_ratio));
				lblEnrichmentSignificance.setText("Enrichment significance: "+String.format("%.3f", significance));
				lblResultsSavedIn.setText("Results saved in: "+work_directory);
			}
		});
		Run_Analysis.setBounds(417, 461, 117, 29);
		tabRunAnalysis.add(Run_Analysis);
		
		
		
		
	}
	
	public boolean ExpressionLoader() throws IOException
	{
		boolean target_loaded = false;
		boolean experiment_loaded = false;
		
		this.experiment_expression_file_loader = new ExpressionFileLoader(experiment_file_path, experiment_sample_identifier_1, experiment_sample_identifier_2);
		experiment_loaded = this.experiment_expression_file_loader.LoadFile();
		
		
		if (main_frame.is_gene_list)
		{
			this.target_gene_list = new GeneListFileLoader(this.target_file_path);
			target_loaded = this.target_gene_list.LoadFile();
		}
		else
		{
			this.target_sample_identifier_1 = target_identifier_1.getText();
			this.target_sample_identifier_2 = target_identifier_2.getText();		
			this.target_expression_file_loader = new ExpressionFileLoader(target_file_path, target_sample_identifier_1, target_sample_identifier_2);
			target_loaded = this.target_expression_file_loader.LoadFile();
		}
		
		if(target_loaded && experiment_loaded)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
	public void CalculateExpression()
	{
		this.experiment_gene_expression = new GeneExpression(experiment_expression_file_loader);
		if (!this.is_gene_list)
		{
			// target input is expression file
			this.target_gene_expression = new GeneExpression(this.target_expression_file_loader);
		}
	}
	
	public void Compare() throws IOException
	{
		
		if (!this.is_gene_list)
		{
			// target input is expression file
			this.comparator = new Comparator(this.experiment_gene_expression, this.target_gene_expression);				   
			this.max_ratio = this.comparator.getMaxRatio();
			this.significance = this.comparator.getSignificance();			
		}
		else
		{
			// target input is gene list
			this.comparator = new Comparator(this.experiment_gene_expression, this.target_gene_list);		
			this.max_ratio = this.comparator.getMaxRatio();
			this.significance = this.comparator.getSignificance();
		}
	}
}