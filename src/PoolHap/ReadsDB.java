package PoolHap; 

import PoolHap.Subtype.Similarity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Random;
import java.util.TreeMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ReadsDB {

	public ReadsDB()
	{
	}
	
	public Collection<Read> getReads()
	{
		return Collections.unmodifiableList(reads);
	}
	
	public Collection<Subtype> getGoldens()
	{
		return Collections.unmodifiableCollection(
			goldenSubtypes.values());
	}

	public int countGolden(String golden)
	{
		int ret = 0;
		for (Read r: reads)
			if (golden.equals(r.getGolden()))
				++ret;
		return ret;
	}
	
	public int countSubtype(Subtype st)
	{
		int ret = 0;
		for (Read r: reads)
			if (st.isCompatible(r))
				++ret;
		return ret;
	}

	public int countSubtypeExclusive(
		Subtype st, Collection<Subtype> subtypes)
	{
		int ret = 0;
		for (Read r: reads)
		{
			if (!st.isCompatible(r))
				continue;
			
			boolean other = false;
			for (Subtype stx: subtypes)
			{
				if (stx == st)
					continue;
				if (stx.isCompatible(r))
				{
					other = true;
					break;
				}
			}
			
			if (!other)
				++ret;
		}
		return ret;
	}
	
	public void load(File file, double sample)
		throws Exception
	{
		reads.clear();
		goldenCounts.clear();
		goldenSubtypes.clear();
		
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
				
				Read r = parseLine(s.trim());
				if (r.equals(null)) continue;
				
				if (rand.nextDouble() > sample)
					continue;
				
				reads.add(r);
				
				String golden = r.getGolden();
				if (golden != null)
				{
					Integer c = goldenCounts.get(golden);
					if (c == null)
						goldenCounts.put(golden, 1);
					else
						goldenCounts.put(golden, c+1);
				}
			}
		}
		catch (Exception e)
		{
			logger.error("{}:{}: {}",
				file, lineNo, e.toString());
			throw e;
		}

		System.out.printf(
			"@name %s, sample %.3f%n",
			file.getName(), sample);
		System.out.printf(
			"@reads %d, lines %d, goldens %d%n",
			reads.size(), lineNo-1, goldenCounts.size());
		
		logger.info(
			"{}: load {} reads from {} lines at rate {}",
			file, reads.size(), lineNo-1, sample);	
		
		goldenCounts.keySet().forEach(g -> {
			goldenSubtypes.put(g, createGoldenSubtype(g));
		});
	}
	
	private static class Rate
	{
		public final int i;
		public final int level;
		public final String levelStr;
		
		public Rate(int i, int level)
		{
			this.i = i;
			this.level = level;
			this.levelStr = levels[level];
		}
		
		public static final int SAME = 0;
		public static final int COMP = 1;
		public static final int ALMO = 2;
		public static final int SIMI = 3;
		public static final int COUL = 4;

		private static final String[] levels = {
			"same as",
			"compatible with",
			"almost same as",
			"similar to",
			"could be"};
	}
	
	public void showSubtypes(String name,
		Collection<Subtype> subtypes, boolean showCounts)
	{
		logger.info(
			"There are {} {} types:",
			subtypes.size(), name);
		
		ArrayList<Subtype> goldens
			= new ArrayList<>(goldenSubtypes.values());
		
		int[] levelCounts = new int[5];
		
		for (Subtype st: subtypes)
		{
			double total = reads.size();
			int g = countGolden(st.getName());
			int cp = countSubtype(st);
			int ex = countSubtypeExclusive(st, subtypes);

			Similarity[] s = Subtype.cmp(st, goldens);
			
			logger.info("  {}: n {}, {}",
				st.getName(), st.size(), Arrays.toString(s));

			if (showCounts)
			{
				Rate rate = rateSubtype(st, s);
				
				if (rate != null)
				{
					logger.info("    {} {} {}",
						rate.levelStr, rate.i,
						goldens.get(rate.i).getName());
					
					++levelCounts[rate.level];
				}
			}

			logger.info(String.format(
				"    freq (%.3f, %.3f, %.3f) (%d, %d, %d)",
				ex/total, g/total, cp/total, ex, g, cp));
		}
		
		if (showCounts)
		{
			logger.info(
				"  Found {} out of {}.",
				Arrays.toString(levelCounts),
				goldens.size());
			System.out.printf(
				"@c_same %d, c_comp %d, c_almo %d, c_simi %d, c_coul %d%n",
				levelCounts[Rate.SAME],
				levelCounts[Rate.COMP],
				levelCounts[Rate.ALMO],
				levelCounts[Rate.SIMI],
				levelCounts[Rate.COUL]);
		}
	}
	
	private Rate rateSubtype(Subtype st, Similarity[] s)
	{
		for (int i = 0; i < s.length; ++i)
		{
			if ((s[i].diff <= 0.01) && (s[i].cover >= 0.99))
			{
				return new Rate(i, Rate.SAME);
			}
		}
		
		for (int i = 0; i < s.length; ++i)
		{
			if (s[i].diff <= 0.01)
			{
				return new Rate(i, Rate.COMP);
			}
		}

		for (int i = 0; i < s.length; ++i)
		{
			if ((s[i].diff <= 0.05) && (s[i].cover >= 0.95))
			{
				return new Rate(i, Rate.ALMO);
			}
		}
		
		for (int i = 0; i < s.length; ++i)
		{
			if ((s[i].diff <= 0.10) && (s[i].cover >= 0.90))
			{
				return new Rate(i, Rate.SIMI);
			}
		}
		
		for (int i = 0; i < s.length; ++i)
		{
			if ((s[i].diff <= 0.15) && (s[i].cover >= 0.85))
			{
				return new Rate(i, Rate.COUL);
			}
		}
		
		return null;
	}
	
	private Subtype createGoldenSubtype(String golden)
	{
		ArrayList<Read> subReads = new ArrayList<>();
		
		for (Read r: reads)
			if (golden.equals(r.getGolden()))
				subReads.add(r);
		
		return new Subtype(golden, subReads);
	}
	
	private Read parseLine(String s)
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
			if (locMut.length != 2) return null; 
				// throw new Exception("invalid location/mutation");
			
			if (pairs.put(Integer.parseInt(locMut[0]),
				Integer.parseInt(locMut[1])) != null)
				throw new Exception("duplicated location");
		}
		assert pairs.size() == locMuts.length;
		
		String[] goldenLoc = (cols.length > 2)?
			cols[2].split(";"): null;
		if ((goldenLoc != null)
			&& (goldenLoc.length < 1))
			throw new Exception("invalid golden column");
		
		return new Read(name, pairs,
			(goldenLoc == null)? null: goldenLoc[0]);
	}
	
	public double totalReads() {
		return reads.size(); 
	}
	
	private final ArrayList<Read> reads = new ArrayList<>();
	private final TreeMap<String, Integer>
		goldenCounts = new TreeMap<>();
	private final TreeMap<String, Subtype>
		goldenSubtypes = new TreeMap<>();
	private final Random rand = new Random(0);
	
	private static final Logger logger
		= LoggerFactory.getLogger(ReadsDB.class);

}
