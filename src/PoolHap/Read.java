package PoolHap; 

import java.util.Map;
import java.util.TreeMap;

public class Read {

	public Read(String name,
			TreeMap<Integer, Integer> pairs,
			String golden)
		{
			this.name = name;
			this.locations = new int[pairs.size()];
			this.mutations = new int[pairs.size()];
			this.golden = golden;
			
			int i = 0;
			for (Map.Entry<Integer, Integer> pair: pairs.entrySet())
			{
				this.locations[i] = pair.getKey();
				this.mutations[i] = pair.getValue();
				++i;
			}
		}
		
		public String getName()
		{
			return name;
		}

		public String getGolden()
		{
			return golden;
		}
		
		public int size()
		{
			return locations.length;
		}
		
		public int getLocation(int i)
		{
			return locations[i];
		}

		public int getMutation(int i)
		{
			return mutations[i];
		}
		
		public int begin()
		{
			return locations[0];
		}
		
		public int end()
		{
			return locations[locations.length-1]+1;
		}
		
		public static boolean isCompatible(Read a, Read b)
		{
			// short cut
			if ((a.begin() >= b.end())
				|| (b.begin() >= a.end()))
				return true;
			
			for (int i = 0, j = 0;
				(i < a.locations.length) && (j < b.locations.length);
				)
			{
				if (a.locations[i] < b.locations[j])
				{
					++i;
				}
				else if (a.locations[i] > b.locations[j])
				{
					++j;
				}
				else // a.locations[i] == b.locations[j]
				{
					if (a.mutations[i] != b.mutations[j])
						return false;
					
					++i;
					++j;
				}
			}
			
			return true;
		}

		private final String name;
		private final int[] locations;
		private final int[] mutations;
		private final String golden;

}
