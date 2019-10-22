// Read gene expression files
import java.io.*;
import java.util.*;


public class ExpressionFileLoader {	
	
	private ArrayList<String> gene_list;
	private String file_path;
	private HashMap sample_dict;
	private String sample_id_1;
	private String sample_id_2;
	private ArrayList<String> sample_name_list_1;
	private ArrayList<String> sample_name_list_2;
	private int num_sample;
	private HashMap<String, HashMap<String, Double>> gene_expression_map;
	
	
	public ExpressionFileLoader(String file_path, String sample_identifier_1, String sample_identifier_2) throws IOException
	{
		// absolute path to the file
		this.file_path = file_path;
		
		// Control identifier
		this.sample_id_1 = sample_identifier_1;
		
		// Treatment identifier
		this.sample_id_2 = sample_identifier_2;
		
		// initialize header list
		this.sample_name_list_1 = new ArrayList<String>();
		this.sample_name_list_2 = new ArrayList<String>();
		
		// initialize gene list
		this.gene_list = null;
		
		// initialize sample dictionary
		this.sample_dict = null;
		
		// <gene_name, <sample_name, expression>>
		this.gene_expression_map = new HashMap<String, HashMap<String, Double>>();
	}
	
	public HashMap<String, HashMap<String, Double>> getGeneExpressionMap() {
		return gene_expression_map;
	}

	public boolean LoadFile()
	{
		BufferedReader buf = null;
		try 
		{
			File file = new File(file_path); 
			buf = new BufferedReader(new FileReader(file));
			
			// read the header line
			String header_line = buf.readLine();
			String[] header_split = header_line.split("\t");
			ArrayList<String> header_list = new ArrayList<String> (Arrays.asList(header_split));
			header_list.remove(0);
			
			// first header is symbol name, sample_1, smaple_2...
			this.num_sample = header_list.size();
			
			this.gene_list = new ArrayList<String>();
			
			this.sample_dict = new HashMap();			
			
			List<List<Double>> temp_all_sample_expression = new ArrayList<List<Double>>();
			
			for (int index=0; index<num_sample; index++)
			{
				temp_all_sample_expression.add(new ArrayList<Double>());
			}
			
			
			String newLine = buf.readLine();
			while(true)
			{
				if(newLine == null)
	            {
	            	// reach EOF
	                break; 
	            }
	            else
	            {
	                String[] expression_line = newLine.split("\t");
	                
	                // ignore all NA genes and repetitive genes
	                if (expression_line[0].equals("NA") || this.gene_list.contains(expression_line[0]))
	                {
	                	newLine = buf.readLine();
	                	continue;
	                }
	                
	                // add gene symbol
	                gene_list.add(expression_line[0]);
	                this.gene_expression_map.put(expression_line[0], new HashMap<String, Double>());
	                
	                // add gene expression
	                for (int col_index=1; col_index<expression_line.length; col_index++)
	                {
	                	temp_all_sample_expression.get(col_index-1).add(Double.valueOf(expression_line[col_index])+1);
	                	this.gene_expression_map.get(expression_line[0]).put(header_list.get(col_index-1), Double.valueOf(expression_line[col_index])+1);
	                }
	            }
				newLine = buf.readLine();
	        }
			
			List<List<Double>> all_sample_expression_1 = new ArrayList<List<Double>>();
			List<List<Double>> all_sample_expression_2 = new ArrayList<List<Double>>();
					
			for (int index=0; index<num_sample; index++)
			{
				if (header_list.get(index).startsWith(this.sample_id_1))
				{
					all_sample_expression_1.add(temp_all_sample_expression.get(index));
					this.sample_name_list_1.add(header_list.get(index));
				}
				else
				{
					all_sample_expression_2.add(temp_all_sample_expression.get(index));
					this.sample_name_list_2.add(header_list.get(index));
				}
			}
			
			this.sample_dict.put(this.sample_id_1, all_sample_expression_1);
			this.sample_dict.put(this.sample_id_2, all_sample_expression_2);
		} catch (Exception e) {
		    e.printStackTrace();
		    return false;
		} finally {
		    try {
		    	buf.close();
		    } catch (IOException e) {
		        e.printStackTrace();
		        return false;
		    }
		}
		return true;
	}
	
	public ArrayList<String> getGeneList()
	{
		return this.gene_list;
	}
	
	public Object getExpr1()
	{
		return this.sample_dict.get(this.sample_id_1);
	}
	
	public Object getExpr2()
	{
		return this.sample_dict.get(this.sample_id_2);
	}
	
	public int getSampleNum()
	{
		return this.num_sample;
	}
	
	public ArrayList<String> getSampleNameList1()
	{
		return this.sample_name_list_1;
	}
	
	public ArrayList<String> getSampleNameList2()
	{
		return this.sample_name_list_2;
	}
	
	public HashMap getAllExpression()
	{
		return this.sample_dict;
	}
	
	
}
