package PoolHap; 

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GoldenDB {

	public GoldenDB()
	{
	}
	
	public Map<String, Subtype> getGoldens()
	{
		return Collections.unmodifiableMap(goldens);
	}

	public void load(File file)
		throws Exception
	{
		goldens.clear();
		
		int lineNo = 1;
		try(
			FileInputStream f = new FileInputStream(file);
			InputStreamReader isr = new InputStreamReader(f, "UTF-8");
			BufferedReader br = new BufferedReader(isr))
		{
			for (;; ++lineNo)
			{
				String s = br.readLine();
				if (s == null)
					break;
				
				Subtype golden = parseLine(s.trim());
				
				goldens.put(golden.getName(), golden);
			}
		}
		catch (Exception e)
		{
			logger.error("{}:{}: {}",
				file, lineNo, e.toString());
			throw e;
		}

		System.out.printf(
			"@gdb_name %s, gdb_lines %d, gdb_size %d%n",
			file.getName(), lineNo-1, goldens.size());
		
		logger.info(
			"{}: load {} goldens from {} lines",
			file, goldens.size(), lineNo-1);	
	}
	
	private Subtype parseLine(String s)
		throws Exception
	{
		String[] cols = s.split("(\\s|:)+");
		if (cols.length < 2)
			throw new Exception("not enough columns");
		
		String name = cols[0];
		
		String[] locMuts = cols[1].trim().split(";");
		if (locMuts.length < 1)
			throw new Exception("no location/mutation");
		
		TreeMap<Integer, Integer> pairs = new TreeMap<>();
		for (String locMutStr: locMuts)
		{
			String[] locMut = locMutStr.split("=");
			if (locMut.length != 2)
				throw new Exception("invalid location/mutation");
			
			if (pairs.put(Integer.parseInt(locMut[0]),
				Integer.parseInt(locMut[1])) != null)
				throw new Exception("duplicated location");
		}
		assert pairs.size() == locMuts.length;
		
		return new Subtype(name, pairs);
	}
	
	private final TreeMap<String, Subtype> goldens = new TreeMap<>();
	
	private static final Logger logger
		= LoggerFactory.getLogger(GoldenDB.class);

}
