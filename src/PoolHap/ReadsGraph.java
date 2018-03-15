package PoolHap; 

import PoolHap.GraphUV.SubGraph;

import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

class ReadsGraph {

	public ReadsGraph(ReadsDB db, ExecutorService workers)
			throws Exception
		{
			Timer t = new Timer();

			this.reads = new ArrayList<>(db.getReads());
			this.reads.sort((l, r) -> l.begin()-r.begin());
			
			this.workers = workers;
			
			this.G = createGraph();

			System.out.printf("@V %d, E %d%n",
				this.G.numVertices(), this.G.numEdges());
			
			logger.info(
				"ReadsGraph: V/E {}/{}, time {}",
				this.G.numVertices(), this.G.numEdges(),
				t.now());
		}
		
		public Read getRead(int v)
		{
			return reads.get(v);
		}
		
		public GraphUV getGraph()
		{
			return G;
		}
		
		public SubGraph getSubGraph(Coloring active)
		{
			IntVector vertices = new IntVector();
			
			vertices.reserve(G.numVertices());
			for (int v = 0; v < G.numVertices(); ++v)
				if (active.isColored(v))
					vertices.pushBack(v);
			
			return G.project(vertices.shrinkToFit());
		}
		
		public ArrayList<Subtype> createSubtypes(int[] colors)
		{
			assert colors.length == G.numVertices();

			ArrayList<Subtype> ret = new ArrayList<>();
			for (int processed = 0; processed < G.numVertices();)
			{
				int c = ret.size();

				ArrayList<Read> subReads = new ArrayList<>();
				for (int v = 0; v < G.numVertices(); ++v)
					if (colors[v] == c)
						subReads.add(reads.get(v));
				
				processed += subReads.size();

				ret.add(new Subtype(
					String.format("Color-%d", c),
					subReads));
			}
			
			return ret;
		}
		
		private GraphUV createGraph()
			throws Exception
		{
			IntVector edges = new IntVector();
			
			ArrayList<Future<int[]>> futures = new ArrayList<>();
			
			final int block = 1000;
			for (int i = 0; i < reads.size(); i += block)
			{
				int begin = i;
				int end = Math.min(i+block, reads.size());
				
				futures.add(workers.submit(
					() -> createEdges(begin, end)));
			}
			
			for (int i = 0; i < futures.size(); ++i)
				edges.append(futures.get(i).get());
			
			return new GraphUV(reads.size(), edges.shrinkToFit());
		}
		
		private int[] createEdges(int begin, int end)
		{
			IntVector edges = new IntVector();
			
			int checked = 0, found = 0;
			for (int from = begin; from < end; ++from)
			{
				Read fromR = reads.get(from);
				for (int to = from+1; to < reads.size(); ++to)
				{
					Read toR = reads.get(to);
					
					// cut off
					if (fromR.end() <= toR.begin())
						break;
					
					// check
					++checked;
					if (!Read.isCompatible(fromR, toR))
					{
						assert !fromR.getGolden().equals(toR.getGolden());
						edges.pushBack(from);
						edges.pushBack(to);
						++found;
					}
				}
			}
			
			logger.debug(
				"ReadsGraph [{}, {}): checked {}, found {}",
				begin, end, checked, found);
			
			return edges.shrinkToFit();
		}
		
		private final ExecutorService workers;
		
		private final ArrayList<Read> reads;
		private final GraphUV G;
		
		private static final Logger logger
			= LoggerFactory.getLogger(ReadsGraph.class);

}
