package PoolHap; 

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

public class Subtype {
	public Subtype(String name, Collection<Read> reads)
	{
		assert name != null;
		this.name = name;
		reads.forEach(r -> add(r));
	}
	
	public Subtype(String name, Map<Integer, Integer> pairs)
	{
		assert name != null;
		this.name = name;
		this.pairs.putAll(pairs);
	}
	
	public String getName()
	{
		return name;
	}
	
	public int size()
	{
		return pairs.size();
	}
	
	public boolean isCompatible(Read r)
	{
		for (int i = 0; i < r.size(); ++i)
		{
			int loc = r.getLocation(i);
			int mut = r.getMutation(i);
			Integer myMut = pairs.get(loc);
			if ((myMut != null) && (myMut != mut))
				return false;
		}
		return true;
	}
	
	@Override
	public String toString()
	{
		StringBuilder builder = new StringBuilder();
		
		builder.append(name).append(": ");
		
		pairs.forEach((l, m) -> {
			builder.append(l).append('=').append(m).append(';');
		});
		
		return builder.toString();
	}
	
	public static void saveSubtypes(	// This is the subroutine that reports the haplotypes. For rjMCMC, we only need the 0/1 in an array.
		File file, Collection<Subtype> subtypes)
		throws Exception
	{
		try(
			FileOutputStream f = new FileOutputStream(file);
			OutputStreamWriter osr = new OutputStreamWriter(f, "UTF-8");
			PrintWriter pw = new PrintWriter(osr))
		{
			for (Subtype subtype: subtypes)
				pw.println(subtype);
		}
	}

	public static HashSet<String> reportSubtypes(Collection<Subtype> initHapList, int[] posArray) throws Exception
	{
		HashSet<String> vefHaplos = new HashSet<String>();
		for (Subtype initHap: initHapList) {
			StringBuilder varComp = new StringBuilder(0);
			for (int p : posArray) {
				if (initHap.pairs.containsKey(p)) varComp.append(initHap.pairs.get(p)); 
				else varComp.append("0"); 
			}
			vefHaplos.add(varComp.toString());
			System.out.println(varComp.toString());
		}
		return vefHaplos; 
	}

	public static class Similarity
	{
		double diff;
		double cover;
		
		public Similarity(double diff, double cover)
		{
			this.diff = diff;
			this.cover = cover;
		}

		@Override
		public String toString()
		{
			return String.format("%.0f/%.0f", diff*100, cover*100);
		}
	}

	public static Similarity cmp(Subtype a, Subtype b)
	{
		int diff = 0, cover = 0;
		for (Map.Entry<Integer, Integer> e: a.pairs.entrySet())
		{
			Integer aMut = a.pairs.get(e.getKey());
			if (aMut == null)
				continue;
			
			++cover;
			
			if (aMut != e.getValue())
				++diff;
		}
		
		return new Similarity(
			diff/(double)b.size(),
			cover/(double)b.size());
	}
	
	public static Similarity[] cmp(Subtype a, Collection<Subtype> bs)
	{
		Similarity[] ret = new Similarity[bs.size()];
		
		int i = 0;
		for (Subtype b: bs)
			ret[i++] = cmp(a, b);
		
		return ret;
	}
	
	private void add(Read r)
	{
		for (int i = 0; i < r.size(); ++i)
		{
			int loc = r.getLocation(i);
			int mut = r.getMutation(i);
			Integer oldMut = pairs.put(loc, mut);
			assert (oldMut == null) || (oldMut == mut);
		}
	}
	
	private final String name;

	private final TreeMap<Integer, Integer>	pairs = new TreeMap<>();
}
