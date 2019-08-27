import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

// for target as gene list file

public class GeneListFileLoader 
{
	private ArrayList<String> gene_list;
	private String file_path;
	
	public GeneListFileLoader(String file_path)
	{
		this.gene_list = new ArrayList<String>();
		this.file_path = file_path;
		// LoadFile();
	}
	
	public boolean LoadFile()
	{
		BufferedReader buf = null;
		try 
		{
			File file = new File(this.file_path); 
			buf = new BufferedReader(new FileReader(file));
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
	            	if (newLine.equals("NA") || this.gene_list.contains(newLine))
	                {
	                	newLine = buf.readLine();
	                	continue;
	                }
	            	else
	            	{
	            		this.gene_list.add(newLine);
	            		newLine = buf.readLine();
	            	}
	            }
			}
			
		}
		catch (IOException e) 
		{
		    e.printStackTrace();
		    return false;
		} 
		finally 
		{
		    try 
		    {
		    	buf.close();
		    } 
		    catch (IOException e) 
		    {
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
	
	public int getGeneNum()
	{
		return this.gene_list.size();
	}
}
