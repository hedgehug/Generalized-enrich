// compare two GeneExpression object or
// compare a GeneExpression object with a GeneList object

import java.awt.Color;
import java.awt.Font;
import java.awt.Paint;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;

// import org.apache.commons.lang3.ArrayUtils;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.ui.RectangleEdge;
import org.tc33.jheatchart.HeatChart;

import net.sf.javaml.utils.ArrayUtils;



public class Comparator 
{
	private GeneExpression experiment_expression;
	private GeneExpression target_expression;
	private GeneListFileLoader target_gene_list;
	private String working_directory;
	
	private double max_ratio;
	private double significance;
	private int window_stride;
	private int initial_start;
	private int permutation_num;
	private int experiment_gene_num;
	private int target_gene_num;
	
	private DefaultCategoryDataset all_line_chart_dataset;
	private HistogramDataset permutation_dataset;
	
	private double[] permutation_max_ratio_list;
	
	private HashMap<String, Double> max_ratio_comparison;
	
	
	public Comparator(GeneExpression experiment_expression, Object target) throws IOException
	{
		this.experiment_expression = experiment_expression;
		this.working_directory = main_frame.work_directory;
		this.max_ratio = 0.0;
		this.significance = 0.0;
		this.experiment_gene_num = experiment_expression.getGeneNum();
		this.permutation_num = main_frame.num_permutation;
		this.window_stride = main_frame.window_stride;
	    this.initial_start = main_frame.initial_window_start;
	    this.max_ratio_comparison = new HashMap<String, Double>();
		
		// initialize figure data set
		this.all_line_chart_dataset = new DefaultCategoryDataset();
		
		
		if(target instanceof GeneExpression)
		{
			this.target_expression = (GeneExpression) target;
			this.target_gene_num = this.target_expression.getGeneNum();
			this.target_gene_list = null;
			DualExpressionCompare();
		}
		else
		{
			this.target_gene_list = (GeneListFileLoader) target;
			this.target_gene_num = this.target_gene_list.getGeneNum();
			this.target_expression = null;
			ExpressionListCompare();
		}
		
		// plot all plots in one figure
		plotAllComparison();		
	}
	
	
	/*
	 * Compare a GeneExpression object with a GeneList object
	 */
	
	public void ExpressionListCompare() throws IOException
	{
		ArrayList<String> experiment_ascending_gene_list = (ArrayList<String>) this.experiment_expression.getSortedGeneList().clone();
		ArrayList<String> experiment_descending_gene_list = (ArrayList<String>) this.experiment_expression.getSortedGeneList().clone();
		Collections.reverse(experiment_descending_gene_list);
		
		ArrayList<String> target_gene_list = (ArrayList<String>) this.target_gene_list.getGeneList().clone();
		
		Compare(experiment_ascending_gene_list, target_gene_list, true, "ascending");
		Compare(experiment_descending_gene_list, target_gene_list, true, "descending");
		
		// permutation test
		permutateTest(experiment_ascending_gene_list, target_gene_list);
	}
	
	
	/*
	 * Compare 2 GeneExpression objects
	 */
	
	public void DualExpressionCompare() throws IOException
	{
		ArrayList<String> experiment_ascending_gene_list = (ArrayList<String>) this.experiment_expression.getSortedGeneList().clone();
		ArrayList<String> experiment_descending_gene_list = (ArrayList<String>) this.experiment_expression.getSortedGeneList().clone();
		Collections.reverse(experiment_descending_gene_list);
		
		ArrayList<String> target_ascending_gene_list = (ArrayList<String>) this.target_expression.getSortedGeneList().clone();
		ArrayList<String> target_descending_gene_list = (ArrayList<String>) this.target_expression.getSortedGeneList().clone();
		Collections.reverse(target_descending_gene_list);
		
		Compare(experiment_ascending_gene_list, target_ascending_gene_list, true, "ascending.vs.ascending");
		Compare(experiment_ascending_gene_list, target_descending_gene_list, true, "ascending.vs.descending");
		Compare(experiment_descending_gene_list, target_descending_gene_list, true, "descending.vs.descending");
		Compare(experiment_descending_gene_list, target_ascending_gene_list, true, "descending.vs.ascending");
		
		// permutation test
		permutateTest(experiment_ascending_gene_list, target_ascending_gene_list);
	}
	
	
	/* 
	 * function Compare() takes two (ranked) gene lists, stored in arraylist of string,
	 * and calculate the odds ratio  of intersection by extending the window from left to right,
	 * the initial size of window and stride of the window extension each time are defined by user
	 * plot the result if required
	 */
	
	public double Compare(ArrayList<String> rank_gene_list_1, ArrayList<String> rank_gene_list_2, boolean if_plotted, String identifier) throws IOException
	{
		Set<String> temp_set_1 = new HashSet<String>();
		Set<String> temp_set_2 = new HashSet<String>();
		int start = 0;
		int exp_end = this.initial_start;
		int tar_end = this.initial_start;
		ArrayList<Double> odds_ratio_list = new ArrayList<Double>(); 
		ArrayList<Integer> gene_position_list = new ArrayList<Integer>(); 
		boolean end_flag = false;
		
		// mark the position of max ratio gene 
		double temp_max_ratio = 0.0;
		int temp_max_ratio_position = exp_end;
		int temp_intersection_num = 0;
		
		while (exp_end<this.experiment_gene_num)
		{
			temp_set_1.addAll(rank_gene_list_1.subList(start, exp_end));
			temp_set_2.addAll(rank_gene_list_2.subList(start, tar_end));
			
			temp_set_1.retainAll(temp_set_2);
			
			int intersection_num = temp_set_1.size();
					
			double odds_ratio = (double)intersection_num / ((double)exp_end / (double)this.experiment_gene_num * (double)exp_end);
			
			if (odds_ratio >= temp_max_ratio)
			{
				temp_max_ratio = odds_ratio;
				temp_max_ratio_position = exp_end;
				temp_intersection_num = intersection_num;
			}
			
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
		
		
		if (if_plotted)
		{
			if (temp_max_ratio > this.max_ratio)
			{
				this.max_ratio = temp_max_ratio;
			}
			
			plotSingleLine(gene_position_list, odds_ratio_list, identifier);
			
			plotHeatMap(rank_gene_list_1, rank_gene_list_2, temp_max_ratio_position, temp_intersection_num, identifier);
			
			// just for return type
			return temp_max_ratio;
		}
		else
		{
			return temp_max_ratio;
		}
	}
	
	/*
	 * Plot the heatmap for gene set  
	 */
	
	public void plotHeatMap(ArrayList<String> rank_gene_list_1, ArrayList<String> rank_gene_list_2, int coordinate, int intersection_num, String identifier) throws IOException
	{
		
		ArrayList<String> header_list = new ArrayList<String>();
		
		
		//ArrayList<ArrayList<Double>> heatmap_data = new ArrayList<ArrayList<Double>>();
		HeatChart heatmap = null;
		
		// write to txt file
		BufferedWriter leading_edge_writer = new BufferedWriter(new FileWriter(main_frame.work_directory+"/leading_edge_gene."+identifier+".txt"));
		
		if (this.target_gene_list == null)
		{
			// if target is gene expression, plot the leading edge
			header_list.addAll(this.target_expression.getSample_name_list_1());
			header_list.addAll(this.target_expression.getSample_name_list_2());
			header_list.addAll(this.experiment_expression.getSample_name_list_1());
			header_list.addAll(this.experiment_expression.getSample_name_list_2());

			int target_sample_num = this.target_expression.getSample_name_list_1().size()+this.target_expression.getSample_name_list_2().size();
			int expr_sample_num = this.experiment_expression.getSample_name_list_1().size()+this.experiment_expression.getSample_name_list_2().size();
			
			int total_sample_num = target_sample_num+expr_sample_num;
			
			double[][] heatmap_data = new double[intersection_num][total_sample_num];
			String[] intersection_gene_list = new String[intersection_num];
			//String[] header_list = new String[total_sample_num]		
			
			ArrayList<String> leading_gene_list_1 = new ArrayList<String> (rank_gene_list_1.subList(0,  coordinate));
			ArrayList<String> leading_gene_list_2 = new ArrayList<String> (rank_gene_list_2.subList(0,  coordinate));
			
			// write header
			leading_edge_writer.write("Gene Name\t"+String.join("\t", header_list)+"\n");
			
			int index = 0;
			int count = 0;
			double temp_min = 0.0;
			double temp_max = 0.0;
			ArrayList<Double> temp_1 = new ArrayList<Double>();
			ArrayList<Double> temp_2 = new ArrayList<Double>();
			while (index<coordinate)
			{
				String gene = leading_gene_list_1.get(index);
				
				if (leading_gene_list_2.contains(gene))
				{
					// write gene name
					leading_edge_writer.write(gene+"\t");					
					
					for (int col_index=0; col_index<target_sample_num; col_index++)
					{
						temp_1.add(this.target_expression.getGeneExpressionMap().get(gene).get(header_list.get(col_index)));
						intersection_gene_list[count] = gene;						
					}
					
					// min-max normalization
					temp_min = Collections.min(temp_1);
					temp_max = Collections.max(temp_1);
					String[] temp_string_1 = new String[target_sample_num];
					for(int i=0; i<temp_1.size(); i++)
					{
						temp_1.set(i, (temp_1.get(i) - temp_min) / (temp_max - temp_min));
					    temp_string_1[i] = Double.toString(temp_1.get(i));
					    heatmap_data[count][i] = temp_1.get(i);
					}						
						
					leading_edge_writer.write(String.join("\t", temp_string_1));
					
					for (int col_index=target_sample_num; col_index<total_sample_num; col_index++)
					{
						temp_2.add(this.experiment_expression.getGeneExpressionMap().get(gene).get(header_list.get(col_index)));
						intersection_gene_list[count] = gene;
					}
					//ArrayUtils.normalize(temp_2);
					// min-max normalization
					temp_min = Collections.min(temp_2);
					temp_max = Collections.max(temp_2);
					String[] temp_string_2 = new String[expr_sample_num];
					for(int i=0; i<temp_2.size(); i++)
					{
						temp_2.set(i, (temp_2.get(i) - temp_min) / (temp_max - temp_min));
					    temp_string_2[i] = String.valueOf(temp_2.get(i));
					    heatmap_data[count][i+temp_1.size()] = temp_2.get(i);
					}						
						
					leading_edge_writer.write(String.join("\t", temp_string_2)+'\n');
					
					temp_1 = new ArrayList<Double>();
					temp_2 = new ArrayList<Double>();

					count++;
				}				
				index++;
			}
			heatmap = new HeatChart(heatmap_data);
			
			heatmap.setHighValueColour(Color.RED);
			heatmap.setLowValueColour(Color.BLUE);
			heatmap.setChartMargin(100);
			
			heatmap.setTitle("Blue-Pink O' Gram");
			heatmap.setYValues(intersection_gene_list);
			heatmap.setXValues(header_list.toArray());
			heatmap.saveToFile(new File(main_frame.work_directory+"/heatmap-"+identifier+".png"));
		}
		else
		{
			// if target is gene list, plot all genes in the gene list
			header_list.addAll(this.experiment_expression.getSample_name_list_1());
			header_list.addAll(this.experiment_expression.getSample_name_list_2());
			
			int total_sample_num = this.experiment_expression.getSample_name_list_1().size()+this.experiment_expression.getSample_name_list_2().size();
			
			double[][] heatmap_data = new double[this.target_gene_list.getGeneList().size()][total_sample_num];
			String[] intersection_gene_list = new String[this.target_gene_list.getGeneList().size()];
			int row_index = 0;
			for (String gene: this.target_gene_list.getGeneList())
			{
				if (rank_gene_list_1.contains(gene))
				{
					//double[] temp_array = new double [total_sample_num];
					for (int col_index=0; col_index<total_sample_num; col_index++)
					{
						heatmap_data[row_index][col_index] = this.experiment_expression.getGeneExpressionMap().get(gene).get(header_list.get(col_index));
						// temp_array[col_index] = this.experiment_expression.getGeneExpressionMap().get(gene).get(header_list.get(col_index));
					}
					ArrayUtils.normalize(heatmap_data[row_index]);
					intersection_gene_list[row_index] = gene;
					row_index++;
				}
			}
			
            heatmap = new HeatChart(heatmap_data);
			
			heatmap.setHighValueColour(Color.RED);
			heatmap.setLowValueColour(Color.BLUE);
			heatmap.setChartMargin(100);
			
			heatmap.setTitle("Blue-Pink O' Gram");
			heatmap.setYValues(intersection_gene_list);
			heatmap.setXValues(header_list.toArray());
			heatmap.saveToFile(new File(main_frame.work_directory+"/heatmap-"+identifier+".png"));
		}	
		leading_edge_writer.close();
	}
	
	
	
	/*
	 * Plot a line plot for any arraylist passed and saved in PNG format
	 * mainly plot max ratio against gene position coordinate
	 */
	
	public void plotSingleLine(ArrayList<Integer> x_coordinate, ArrayList<Double> y_coordinate, String identifier) throws IOException
	{
		double temp_max_ratio = Collections.max(y_coordinate);
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
		
		// TODO refine the figure quality, maybe change dpi
		
		// change the font
		CategoryPlot plot = lineChartObject.getCategoryPlot();
		CategoryAxis x_axis = plot.getDomainAxis();
		x_axis.setLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		x_axis.setTickLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		plot.getRangeAxis().setLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		plot.getRangeAxis().setTickLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		
		// change ticks
		NumberAxis y_axis = (NumberAxis)plot.getRangeAxis(); 
		y_axis.setTickUnit(new NumberTickUnit(1)); 
		
		int width = 1000;    /* Width of the image */
	    int height = 800;   /* Height of the image */ 
	    lineChartObject.removeLegend();
	    TextTitle legendText = new TextTitle("Max Ratio = "+Double.toString(temp_max_ratio));
	    legendText.setPosition(RectangleEdge.BOTTOM);
	    lineChartObject.addSubtitle(legendText);
	    File lineChart = new File(main_frame.work_directory+"/"+identifier+".png"); 
	    ChartUtilities.saveChartAsPNG(lineChart ,lineChartObject, width ,height);
	    
	    this.max_ratio_comparison.put(identifier, temp_max_ratio);
	}
	
	
	public void plotAllComparison() throws IOException
	{
		JFreeChart all_lineChartObject = ChartFactory.createLineChart(
				 "Enrichment", "Gene Number", "Odds Ratio",
				 this.all_line_chart_dataset,PlotOrientation.VERTICAL,
		         true,true,false);
		int width = 1000;    /* Width of the image */
	    int height = 800;   /* Height of the image */ 
	    
	    // change font, ticks
	    CategoryPlot plot = all_lineChartObject.getCategoryPlot();
		CategoryAxis x_axis = plot.getDomainAxis();
		x_axis.setLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		x_axis.setTickLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		plot.getRangeAxis().setLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		plot.getRangeAxis().setTickLabelFont(new Font("TimesRoman", Font.BOLD, 25));
		NumberAxis y_axis = (NumberAxis)plot.getRangeAxis(); 
		y_axis.setTickUnit(new NumberTickUnit(1));   
	    
	    
	    // add customized legend
	    //all_lineChartObject.removeLegend();
	    String temp_legned = "Max ratio = "+Double.toString(this.max_ratio);
	    TextTitle legendText = new TextTitle(temp_legned);
	    legendText.setPosition(RectangleEdge.BOTTOM);
	    all_lineChartObject.addSubtitle(legendText);
	    temp_legned = "p-value = "+Double.toString(this.significance);
	    legendText = new TextTitle(temp_legned);
	    legendText.setPosition(RectangleEdge.BOTTOM);
	    all_lineChartObject.addSubtitle(legendText);	    
	    
	    
	    File lineChart = new File(main_frame.work_directory+"/all_result.png");
	    ChartUtilities.saveChartAsPNG(lineChart , all_lineChartObject, width ,height);
	}
	
	
	/*
	 * Conduct permutation test and calculate significances 
	 */
	
	public void permutateTest(ArrayList<String> rank_gene_list_1, ArrayList<String> rank_gene_list_2) throws IOException
	{
		this.permutation_dataset = new HistogramDataset();
        permutation_dataset.setType(HistogramType.RELATIVE_FREQUENCY);
		
		ArrayList<String> temp_gene_list_1 = (ArrayList<String>) rank_gene_list_1.clone();
		ArrayList<String> temp_gene_list_2 = (ArrayList<String>) rank_gene_list_2.clone();
		
        this.permutation_max_ratio_list = new double[this.permutation_num];
		
		String permutation_result_file_name = main_frame.work_directory+"/permutation_result-"+Integer.toString(this.permutation_num)+".txt";
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
				this.permutation_max_ratio_list[round] = Compare(temp_gene_list_1, temp_gene_list_2, false, Integer.toString(round));
			}
		}
		
        this.permutation_dataset.addSeries("Histogram", this.permutation_max_ratio_list, 20);
		
		JFreeChart permutation_histogram = ChartFactory.createHistogram("Permutation max odds ratio distribution", "number", 
				"Max odds ratio", permutation_dataset, PlotOrientation.VERTICAL, true,true,false);
		
		Paint[] paintArray = null;
	    paintArray = new Paint[1];
        paintArray[0] = new Color(0x800000ff, true);
				
		permutation_histogram.removeLegend();
		int width = 1000;    /* Width of the image */
	    int height = 800;   /* Height of the image */ 
	    ValueMarker marker = new ValueMarker(this.max_ratio);
	    marker.setPaint(Color.RED);
	    XYPlot plot = (XYPlot) permutation_histogram.getPlot();
	    
	    // set histogram bars color to blue
	    plot.setDrawingSupplier(new DefaultDrawingSupplier(paintArray,
	            DefaultDrawingSupplier.DEFAULT_PAINT_SEQUENCE,
	            DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE,
	            DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE,
	            DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE));
	    
	    plot.addDomainMarker(marker);
	    
	    File hist = new File(main_frame.work_directory+"/permutation_result.png"); 
	    ChartUtilities.saveChartAsPNG(hist , permutation_histogram, width ,height);
	    
	    Arrays.sort(this.permutation_max_ratio_list);
	    
	    BufferedWriter writer = new BufferedWriter(new FileWriter(permutation_result_file_name));
	    for (int index=0; index<this.permutation_max_ratio_list.length; index++)
	    {
	    	writer.write(Double.toString(this.permutation_max_ratio_list[index])+"\n");
	    }
	    writer.close();
	    
	    int counter = 0;
	    while (counter<this.permutation_max_ratio_list.length && this.permutation_max_ratio_list[counter]<this.max_ratio)
	    {
	    	counter++;
	    }
	    this.significance = 1 - (double) counter / (double) this.permutation_num;
	}
	
	public double getMaxRatio()
	{
		return this.max_ratio;
	}
	
	public double getSignificance()
	{
		return this.significance;
	}
	
	public void writeHTML() throws IOException
	{
		BufferedWriter html_writer = new BufferedWriter(new FileWriter(main_frame.work_directory+"/report.html"));
	    for (int index=0; index<this.permutation_max_ratio_list.length; index++)
	    {
	    	// write HTML header
	    	html_writer.write("<!DOCTYPE html>\n<html>\n<body>\n");
	    	
	    	// TODO write the body of HTML, index.html for reference
	    	html_writer.write("<h1>Result of Generalized Gene Set Enrichment Analysis</h1>");
	    	
	    	// write HTML footer
	    	html_writer.write("</body>\n</html>\n");
	    }
	    html_writer.close();
	}

}
