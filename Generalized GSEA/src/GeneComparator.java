import java.awt.Color;
import java.awt.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

// take two GeneExpression Objects and make comparison


public class GeneComparator 
{
	
	private int window_stride;
	private int initial_start;
	private int experiment_gene_num;
	private int target_gene_num;
	private double max_odds_ratio;
	private DefaultCategoryDataset all_line_chart_dataset;
	private HistogramDataset permutation_dataset;
	private int permutation_num;
	private double[] permutation_max_ratio_list;
	private double significance;
	
	public GeneComparator(GeneExpression experiment_gene_expression, GeneExpression target_gene_expression, int stride, int initial_window_start, int permutation_num) throws IOException
	{
		this.window_stride = stride;
		this.initial_start = initial_window_start;
		this.experiment_gene_num = experiment_gene_expression.getGeneNum();
		this.target_gene_num = target_gene_expression.getGeneNum();
		this.max_odds_ratio = 0;
		this.permutation_num = permutation_num;
		
		ArrayList<String> experiment_ascending_gene_list = (ArrayList<String>) experiment_gene_expression.getSortedGeneList().clone();
		ArrayList<String> experiment_descending_gene_list = (ArrayList<String>) experiment_gene_expression.getSortedGeneList().clone();
		Collections.reverse(experiment_descending_gene_list);
		
		ArrayList<String> target_ascending_gene_list = (ArrayList<String>) target_gene_expression.getSortedGeneList().clone();
		ArrayList<String> target_descending_gene_list = (ArrayList<String>) target_gene_expression.getSortedGeneList().clone();
		Collections.reverse(target_descending_gene_list);

		
		this.all_line_chart_dataset = new DefaultCategoryDataset();
		
		Compare(experiment_ascending_gene_list, target_ascending_gene_list, "ascending.vs.ascending");
		Compare(experiment_ascending_gene_list, target_descending_gene_list, "ascending.vs.descending");
		Compare(experiment_descending_gene_list, target_descending_gene_list, "descending.vs.descending");
		Compare(experiment_descending_gene_list, target_ascending_gene_list, "descending.vs.ascending");
		
		
		
		JFreeChart all_lineChartObject = ChartFactory.createLineChart(
				 "Enrichment", "Gene Number", "Odds Ratio",
				 this.all_line_chart_dataset,PlotOrientation.VERTICAL,
		         true,true,false);
		
		int width = 1000;    /* Width of the image */
	    int height = 800;   /* Height of the image */ 
	    File lineChart = new File(main_frame.work_directory+"/all_result.png"); 
	    ChartUtilities.saveChartAsPNG(lineChart , all_lineChartObject, width ,height);
	    
	    randomTest(experiment_ascending_gene_list, target_ascending_gene_list);
	}
	
	public void Compare(ArrayList<String> rank_gene_list_1, ArrayList<String> rank_gene_list_2, String identifier) throws IOException
	{
		Set<String> temp_set_1 = new HashSet<String>();
		Set<String> temp_set_2 = new HashSet<String>();
		int start = 0;
		int exp_end = this.initial_start;
		int tar_end = this.initial_start;
		ArrayList<Double> odds_ratio_list = new ArrayList<Double>(); 
		ArrayList<Integer> gene_position_list = new ArrayList<Integer>(); 
		boolean end_flag = false;
		
		while (exp_end<this.experiment_gene_num)
		{
			temp_set_1.addAll(rank_gene_list_1.subList(start, exp_end));
			temp_set_2.addAll(rank_gene_list_2.subList(start, tar_end));
			
			temp_set_1.retainAll(temp_set_2);
			
			int intersection_num = temp_set_1.size();
					
			double odds_ratio = (double)intersection_num / ((double)exp_end / (double)this.experiment_gene_num * (double)exp_end);
			
			odds_ratio_list.add(odds_ratio);
			gene_position_list.add(exp_end);
			if (end_flag)
			{
				break;
			}
			exp_end += this.window_stride;
			tar_end += this.window_stride;
			
			if (tar_end>=this.target_gene_num)
			{
				tar_end = this.target_gene_num-1;
			}
			
			if (exp_end>=this.experiment_gene_num)
			{
				exp_end = this.experiment_gene_num-1;
				end_flag = true;			
			}
		}
		
		if (Collections.max(odds_ratio_list) > this.max_odds_ratio)
		{
			this.max_odds_ratio = Collections.max(odds_ratio_list);
		}
		
		plotSeperate(gene_position_list, odds_ratio_list, identifier);
	}
	
	public void plotSeperate(ArrayList<Integer> x_coordinate, ArrayList<Double> y_coordinate, String identifier) throws IOException
	{
		DefaultCategoryDataset line_chart_dataset = new DefaultCategoryDataset();
		int point_num = x_coordinate.size();
		
		for (int index=0; index<point_num; index++)
		{
			line_chart_dataset.addValue(y_coordinate.get(index), identifier, x_coordinate.get(index));
			this.all_line_chart_dataset.addValue(y_coordinate.get(index), identifier, x_coordinate.get(index));;
		}
		
		JFreeChart lineChartObject = ChartFactory.createLineChart(
				 identifier, "Gene Number", "Odds Ratio",
		         line_chart_dataset,PlotOrientation.VERTICAL,
		         true,true,false);
		
		int width = 1000;    /* Width of the image */
	    int height = 800;   /* Height of the image */ 
	    File lineChart = new File(main_frame.work_directory+"/"+identifier+".png"); 
	    ChartUtilities.saveChartAsPNG(lineChart ,lineChartObject, width ,height);
	}
	
	public void randomTest(ArrayList<String> rank_gene_list_1, ArrayList<String> rank_gene_list_2) throws IOException
	{
		this.permutation_dataset = new HistogramDataset();
		permutation_dataset.setType(HistogramType.RELATIVE_FREQUENCY);
		
		ArrayList<String> temp_gene_list_1 = (ArrayList<String>) rank_gene_list_1.clone();
		ArrayList<String> temp_gene_list_2 = (ArrayList<String>) rank_gene_list_2.clone();
		
		this.permutation_max_ratio_list = new double[this.permutation_num];
		
		String permutation_result_file_name = main_frame.work_directory+"/permutation_result.txt";
		File temp_permutation = new File(permutation_result_file_name);
		if (temp_permutation.exists())
		{
			BufferedReader buf = new BufferedReader(new FileReader(temp_permutation));
			String newLine = buf.readLine();
			int row_count = 0;
			while(true)
			{
				if(newLine == null)
	            {
	            	// reach EOF
	                break; 
	            }
	            else
	            {
	            	this.permutation_max_ratio_list[row_count] = Double.parseDouble(newLine);
	                newLine = buf.readLine();
	                row_count++;
	            }
			}
			buf.close();
		}
		else 
		{
			for (int round=0; round<this.permutation_num; round++)
			{
				Collections.shuffle(temp_gene_list_1);
				this.permutation_max_ratio_list[round] = randomCompare(temp_gene_list_1, temp_gene_list_2, Integer.toString(round));
			}
		}
		
		
		permutation_dataset.addSeries("Histogram", this.permutation_max_ratio_list, 20);
		
		JFreeChart permutation_histogram = ChartFactory.createHistogram("Permutation max odds ratio distribution", "number", 
				"Max odds ratio", permutation_dataset, PlotOrientation.VERTICAL, true,true,false);
		
		permutation_histogram.removeLegend();
		int width = 1000;    /* Width of the image */
	    int height = 800;   /* Height of the image */ 
	    ValueMarker marker = new ValueMarker(this.max_odds_ratio);
	    marker.setPaint(Color.RED);
	    XYPlot plot = (XYPlot) permutation_histogram.getPlot();
	    plot.addDomainMarker(marker);
	    File hist = new File(main_frame.work_directory+"/permutation_result.png"); 
	    ChartUtilities.saveChartAsPNG(hist , permutation_histogram, width ,height);
	    
	    // java.util.List<double[]> temp_max_ratio_list = Arrays.asList(this.permutation_max_ratio_list);
	    Arrays.sort(this.permutation_max_ratio_list);
	    
	    BufferedWriter writer = new BufferedWriter(new FileWriter(main_frame.work_directory+"/permutation_result.txt"));
	    for (int index=0; index<this.permutation_max_ratio_list.length; index++)
	    {
	    	writer.write(Double.toString(this.permutation_max_ratio_list[index])+"\n");
	    }
	    writer.close();
	    
	    int counter = 0;
	    while (this.permutation_max_ratio_list[counter]<this.max_odds_ratio)
	    {
	    	counter++;
	    }
	    this.significance = 1 - (double) counter / (double) this.permutation_num;
	}
	
	public Double randomCompare(ArrayList<String> rank_gene_list_1, ArrayList<String> rank_gene_list_2, String identifier) throws IOException
	{
		Set<String> temp_set_1 = new HashSet<String>();
		Set<String> temp_set_2 = new HashSet<String>();
		int start = 0;
		int exp_end = this.initial_start;
		int tar_end = this.initial_start;
		boolean end_flag = false;
		ArrayList<Double> temp_max_ratio = new ArrayList<Double>();
		
		while (exp_end<this.experiment_gene_num)
		{
			temp_set_1.addAll(rank_gene_list_1.subList(start, exp_end));
			temp_set_2.addAll(rank_gene_list_2.subList(start, tar_end));
			
			temp_set_1.retainAll(temp_set_2);
			
			int intersection_num = temp_set_1.size();
					
			double odds_ratio = (double)intersection_num / ((double)exp_end / (double)this.experiment_gene_num * (double)exp_end);
			
			temp_max_ratio.add(odds_ratio);
			// this.random_line_chart_dataset.addValue((double)odds_ratio, identifier, Integer.toString(exp_end));
			
			if (end_flag)
			{
				break;
			}
			
			exp_end += this.window_stride;
			tar_end += this.window_stride;
			
			if (tar_end>=this.target_gene_num)
			{
				tar_end = this.target_gene_num-1;
			}
			
			if (exp_end>=this.experiment_gene_num)
			{
				// ends comparison when sliding to the end of the experiment gene list
				exp_end = this.experiment_gene_num-1;
				end_flag = true;			
			}
		}
		return Collections.max(temp_max_ratio);
	}
	
	public double getSignificance()
	{
		return this.significance;
	}
	
	public double getMaxRatio()
	{
		return this.max_odds_ratio;
	}
}
