package PoolHap; 

import PoolHap.GraphUV.SubGraph;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SubsetReduction {

	public SubsetReduction(GraphUV G)
	{
		this.G = G;
		
		final int numV = G.numVertices();
		final int numE = G.numEdges();
		
		this.active = new Coloring(numV);
		for (int v = 0; v < numV; ++v)
			this.active.setColor(v);
		
		this.begins = new int[numV];
		this.ends = new int[numV];
		this.neighbors = new int[numE*2];

		this.activeV = new IntVector();
		this.activeV.reserve(numV);
		
		this.hints = new AtomicIntegerArray(numV);
		for (int v = 0; v < numV; ++v)
			this.hints.lazySet(v, -2);
	}
	
	private void buildActiveSubgraph()
	{
		activeV.clear();
		for (int v = 0, cur = 0; v < G.numVertices(); ++v)
		{
			if (!active.isColored(v))
			{
				// as precaution
				begins[v] = ends[v] = cur;
				continue;
			}
			
			activeV.pushBack(v);
			
			begins[v] = cur;
			
			for (int iv = G.neighborsBegin(v);
				iv != G.neighborsEnd(v); ++iv)
			{
				int u = G.neighborVertex(iv);
				if (active.isColored(u))
					neighbors[cur++] = u;
			}
			
			ends[v] = cur;
			
			// as needed by subset check
			Arrays.sort(neighbors, begins[v], ends[v]);
		}
	}
	
	private void deactivate(int v, int hint)
	{
		logger.debug("v {}, hint {}", v, hint);
		
		assert active.isColored(v);
		active.resetColor(v);
		
		assert (hint == -1) || active.isColored(hint);
		hints.lazySet(v, hint);
	}
	
	public int[] resolveColors(int[] subColors)
	{
		assert activeV.size() == subColors.length;

		int[] colors = new int[G.numVertices()];
		Arrays.fill(colors, -1);

		// color active vertices first
		for (int subV = 0; subV < subColors.length; ++subV)
		{
			int v = activeV.get(subV);
			logger.debug(
				"v {}, subV {}, color {}",
				v, subV, subColors[subV]);
			colors[v] = subColors[subV];
		}
		
		// then use hints
		for (int v = 0; v < G.numVertices(); ++v)
		{
			// All active vertices should be colored already.
			if (active.isColored(v))
			{
				assert colors[v] != -1;
				continue;
			}
			
			// already colored in DFS
			if (colors[v] != -1)
				continue;
			
			dfsColoring(v, colors);
		}
		
		return colors;
	}
	
	private void dfsColoring(int v, int[] colors)
	{
		int u = hints.get(v);
		
		if (u == -1)
		{
			logger.debug("v {}, hint -1, color 0", v);
			colors[v] = 0;
			return;
		}
		
		if (colors[u] == -1)
		{
			assert !active.isColored(u);
			dfsColoring(u, colors);
		}
		
		logger.debug(
			"v {}, hint {}, color {}",
			v, u, colors[u]);
		
		colors[v] = colors[u];
	}
	
	public SubGraph apply()
	{
		Timer t = new Timer();
		
		for (int round = 1;; ++round)
		{
			totalChecks = totalLoops = 0;
			
			buildActiveSubgraph();
			
			logger.info(
				"Subset {}: V/E {}/{}",
				round, activeV.size(),
				ends[G.numVertices()-1]/2);

			int reduced = oneRound();
			
			logger.info(
				"Subset {}: t {}, reduced {}, {}/{}",
				round, t.now(), reduced,
				totalLoops, totalChecks);
			
			if (reduced == 0)
				break;
		}
		
		// Since the loop exits when there is no change,
		// there is no need to re-calculate active vertices.
		SubGraph sG = G.project(activeV.shrinkToFit());
		
		System.out.printf(
			"@V_ssr %d, E_ssr %d, real_ssr %.3f%n",
			sG.numVertices(), sG.numEdges(), t.now());

		logger.info(
			"Subset done: t {}, V/E {}/{}",
			t.now(), sG.numVertices(), sG.numEdges());
		
		return sG;
	}

	private int oneRound()
	{
		int ret = 0;
		
		// neighbors/distance-2 neighbors
		Coloring nbs = new Coloring(G.numVertices());
		Coloring nbs2 = new Coloring(G.numVertices());
		IntVector nbs2V = new IntVector();
		
		totalChecks = totalLoops = 0;
		for (int i = 0; i < activeV.size(); ++i)
		{
			int v = activeV.at(i);
			
			// already deactivated?
			if (!active.isColored(v))
			{
				continue;
			}
			
			// 0 degree? color as you like
			if (ends[v] == begins[v])
			{
				deactivate(v, -1);
				++ret;
				continue;
			}
			
			// init neighbors
			nbs.resetAll();
			nbs2.resetAll();
			nbs2V.clear();
			for (int iv = begins[v]; iv != ends[v]; ++iv)
			{
				int x = neighbors[iv];
				if (!active.isColored(x))
					continue;
				
				nbs.setColor(x);
				
				for (int ix = begins[x]; ix != ends[x]; ++ix)
				{
					int u = neighbors[ix];
					if ((u == v) // no self check
						|| !active.isColored(u) // not active
						|| nbs2.isColored(u)) // will check
						continue;
					
					nbs2V.pushBack(u);
					nbs2.setColor(u);
				}
			}
			
			// check them
			for (int j = 0; j < nbs2V.size(); ++j)
			{
				int u = nbs2V.at(j);
				
				if (isSubset(u, nbs))
				{
					deactivate(u, v);
					++ret;
				}
			}
		}
		
		return ret;
	}
	
	private boolean isSubset(int u, Coloring nbs)
	{
		++totalChecks;
		
		for (int iu = begins[u]; iu < ends[u]; ++iu)
		{
			++totalLoops;
		
			int x = neighbors[iu];
			
			// ignore
			if (!active.isColored(x))
				continue;
			
			// not v's neighbor
			if (!nbs.isColored(x))
				return false;
		}
		
		return true;
	}
	
	private final GraphUV G;

	// active subgraph
	private final Coloring active;
	private final int[] begins;
	private final int[] ends;
	private final int[] neighbors;
	private final IntVector activeV;

	// hints
	private final AtomicIntegerArray hints;
	
	private long totalChecks, totalLoops;
	
	private static final Logger logger
		= LoggerFactory.getLogger(SubsetReduction.class);

}
