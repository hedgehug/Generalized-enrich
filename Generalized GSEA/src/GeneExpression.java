// gene expression holder

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;


public class GeneExpression {
	
	private HashMap<String, double[]> sample_stats_1;
	private HashMap<String, double[]> sample_stats_2;
	private int gene_num;
	private int sample_num;
	private ArrayList<String> gene_list;
	// sort in ascending order
	private ArrayList<String> sorted_gene_list;
	private ArrayList<Double> sorted_expression_list;
	
	
	public GeneExpression(ExpressionFileLoader file_dataset)
	{
		this.sample_stats_1 = Calculate((List<List<Double>>) file_dataset.getExpr1());		
		this.sample_stats_2 = Calculate((List<List<Double>>) file_dataset.getExpr2());
		this.gene_list = file_dataset.getGeneList();
		this.sorted_gene_list = new ArrayList<String> ();
		this.sorted_expression_list = new ArrayList<Double> ();
		
		Rank_gene();
	}
	
	public HashMap<String, double[]> Calculate(List<List<Double>> sample_expression)
	{		
		this.gene_num = sample_expression.get(0).size();
		this.sample_num = sample_expression.size();
		
		double[] sample_mean = new double[this.gene_num];
		double[] sample_std = new double[this.gene_num];
		// List<Double> sample_std = new ArrayList<Double>();
		
		for (int row_index=0; row_index<this.gene_num; row_index++)
		{
			// calculate mean
			double temp_sum = 0;
			for (int col_index=0; col_index<this.sample_num; col_index++)
			{
				temp_sum += sample_expression.get(col_index).get(row_index);
			}
			double mean = temp_sum / this.sample_num;
			if (mean==0)
			{
				mean = 1;
			}
			sample_mean[row_index] = mean;
			
			// calculate standard deviation
			temp_sum = 0;
			for (int col_index=0; col_index<this.sample_num; col_index++)
			{
				temp_sum += Math.pow((sample_expression.get(col_index).get(row_index)-mean), 2);
			}
			double std = Math.max(Math.sqrt(temp_sum / (this.sample_num-1)), 0.2*Math.abs(mean));
			sample_std[row_index] = std;			
		}
		
		HashMap<String, double[]> stats = new HashMap<String, double[]>();
		stats.put("mean", sample_mean);
		stats.put("std", sample_std);
		
		return stats;
	}
	
	public void Rank_gene()
	{
		double[] temp_expression = new double[this.gene_num];
		for (int row_index=0; row_index<this.gene_num; row_index++)
		{
			temp_expression[row_index] = ((sample_stats_2.get("mean")[row_index]-sample_stats_1.get("mean")[row_index]+0.0001)/
					(sample_stats_2.get("std")[row_index]+sample_stats_1.get("std")[row_index]+0.0001));
		}
		
		HashMap<String, Double> temp_gene_expression_dict = new HashMap<String, Double>();
		
		for (int row_index=0; row_index<this.gene_num; row_index++)
		{
			temp_gene_expression_dict.put(gene_list.get(row_index), temp_expression[row_index]);
		}
		
		sortByValue(temp_gene_expression_dict);
	}
	
	public void sortByValue(HashMap<String, Double> gene_expression) 
    { 
        // Create a list from elements of HashMap 
        List<Map.Entry<String, Double> > list = 
               new LinkedList<Map.Entry<String, Double> >(gene_expression.entrySet()); 
  
        // Sort the list 
        Collections.sort(list, new Comparator<Map.Entry<String, Double> >() { 
            public int compare(Map.Entry<String, Double> o1,  
                               Map.Entry<String, Double> o2) 
            { 
                return (o1.getValue()).compareTo(o2.getValue()); 
            } 
        }); 
        
        for (Map.Entry<String, Double> aa : list) { 
            // temp.put(aa.getKey(), aa.getValue());
        	this.sorted_gene_list.add(aa.getKey());
        	this.sorted_expression_list.add(aa.getValue());
        } 
    }
	
	public ArrayList<String> getSortedGeneList()
	{
		return this.sorted_gene_list;
	}
	
	public ArrayList<Double> getSortedExpressionList()
	{
		return this.sorted_expression_list;
	}
	
	public int getGeneNum()
	{
		return this.gene_num;
	}
}
