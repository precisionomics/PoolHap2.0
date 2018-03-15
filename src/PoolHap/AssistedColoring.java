package PoolHap; 

import PoolHap.Timer;

import java.util.Arrays;
import java.util.TreeSet;
import java.util.function.Function;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AssistedColoring {

	public AssistedColoring(String name,
			GraphUV G, int[] clique, GoldenDB gdb,
			Function<Integer, Read> reads)
		{
			assert clique.length <= G.numVertices();
			assert clique.length > 0;
			
			this.name = name;
			this.G = G;
			this.clique = clique;
			this.gdb = gdb;
			this.reads = reads;
			
			this.exclude = new Coloring(G.numVertices());
			this.mc = new MaxClique4(G);
			
			this.colors = new int[G.numVertices()];
			
			this.nbColors = new Coloring(G.numVertices());
		}
		
		public int[] run()
		{
			collectCandidates();
			
			generateColors();
			
			return colors;
		}
		
		private void generateColors()
		{
			Arrays.fill(colors, -1);
			
			int numColored = 0;
			int c = 0;
			for (GoldenInfo info: candidates)
			{
				for (int v = 0; v < G.numVertices(); ++v)
				{
					if (colors[v] != -1)
						continue;
					
					Read r = reads.apply(v);
					
					if (!info.golden.isCompatible(r))
						continue;
					
					colors[v] = c;
					++numColored;
				}
				
				++c;
			}
			
			logger.info(
				"Assist: {} vertices colored, remaining {}",
				numColored, G.numVertices()-numColored);
			
			if (numColored != G.numVertices())
			{
				for (int v = 0; v < G.numVertices(); ++v)
				{
					if (colors[v] != -1)
						continue;
					
					colorOne(v);

					Read r = reads.apply(v);
					
					logger.debug(
						"v {}, color {}, {} w/ {}",
						v, colors[v], r.getName(), r.getGolden());
				}
			}
		}
		
		private void collectCandidates()
		{
			candidates.clear();
			
			gdb.getGoldens().forEach((name, golden) -> {
				
				exclude.resetAll();
				
				int count = 0;
				for (int v = 0; v < G.numVertices(); ++v)
				{
					Read r = reads.apply(v);
					
					if (!golden.isCompatible(r))
						continue;
					
					++count;
					exclude.setColor(v);
				}
							
				int target = clique.length-1;
				
				logger.debug(
					"Assist: try {}, count {}, target {}",
					name, count, target);

				int[] exc = mc.solve(exclude, target, target+1, false);
				if (exc.length == target)
				{
					candidates.add(new GoldenInfo(golden, count));
					logger.info(
						"Assist: find {}, count {}, t {}",
						name, count, t.now());
				}
				else
				{
					logger.info(
						"Assist: ignore {}, count {}, t {}",
						name, count, t.now());
				}
			});
			
			logger.info(
				"Assist: found {} candidates for {}",
				candidates.size(), name);
		}
		
		private void colorOne(int v)
		{
			nbColors.resetAll();
			
			for (int iv = G.neighborsBegin(v);
				iv != G.neighborsEnd(v); ++iv)
			{
				int u = G.neighborVertex(iv);
				
				// no color, ignore
				if (colors[u] == -1)
					continue;
				
				nbColors.setColor(colors[u]);
			}
			
			int c = candidates.size(); // use new color
			for (; c < nbColors.size(); ++c)
			{
				if (!nbColors.isColored(c))
					break;
			}
			colors[v] = c;
		}
		
		private final String name;
		private final GraphUV G;
		private final int[] clique;
		private final GoldenDB gdb;
		private final Function<Integer, Read> reads;
		
		private final Coloring exclude;
		private final MaxClique4 mc;
		
		private static class GoldenInfo
			implements Comparable<GoldenInfo>
		{
			public final Subtype golden;
			public final int count;

			public GoldenInfo(Subtype golden, int count)
			{
				this.golden = golden;
				this.count = count;
			}

			@Override
			public int compareTo(GoldenInfo r)
			{
				if (count > r.count)
					return -1;
				else if (count < r.count)
					return 1;
				else
					return golden.getName().compareTo(
						r.golden.getName());
			}
		}
		
		private final TreeSet<GoldenInfo> candidates = new TreeSet<>(); 

		private final int[] colors;
		
		// colors used by neighbors
		private final Coloring nbColors;
		
		private final Timer t = new Timer();
		
		private static final Logger logger
			= LoggerFactory.getLogger(AssistedColoring.class);
}
