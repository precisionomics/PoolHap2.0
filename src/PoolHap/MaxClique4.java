package PoolHap; 

import java.util.Arrays;
import java.util.LinkedList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MaxClique4 {

	public MaxClique4(GraphUV G)
	{
		this.G = G;
		this.clique = new IntVector();
		this.best = new IntVector();
		
		this.subVC = new VCPruner4(G);
		
		this.initBuffer = new int[G.numVertices()];
		
		this.buffersArray = new int[G.numEdges()*2];
		this.buffersOffsets = new int[G.numVertices()];

		// System.out.println(G.numVertices());
		// System.out.println(buffersOffsets.length);

		buffersOffsets[0] = 0;
		for (int v = 1; v < G.numVertices(); ++v)
		{
			buffersOffsets[v] = buffersOffsets[v-1]
				+G.numNeighbors(v-1);
		}
		
		this.depths = new int[G.numVertices()];
	}
	
	public int[] solve()
	{
		return solve(null, 0, -1, true);
	}
	
	public int[] solve(Coloring exclude,
		int lb, int ub, boolean verbose)
	{
		if (lb < 1)
			lb = 1;
		
		assert (ub == -1) || (lb <= ub);
		
		init(exclude, lb, verbose);
		
		l1Count = 0;
		for (; !stack.isEmpty();)
		{
			if (best.size() == ub)
				break;
			
			int v = next(lb, verbose);
			if (v != -1)
			{
				push(v, lb, verbose);
			}
			else
			{
				pop();
			}
		}
		
		if (verbose)
			logger.info(
				"Clique Done: t {}, best {}, lb {}, ub {}",
				t.now(), best.size(), lb, ub);
		else
			logger.debug(
				"Clique Done: t {}, best {}, lb {}, ub {}",
				t.now(), best.size(), lb, ub);
		
		if (best.size() < lb)
			return null;
		else
			return Arrays.copyOf(best.data(), best.size());
	}
	
	private void init(Coloring exclude, int lb, boolean verbose)
	{
		clique.clear();
		best.clear();
		stack.clear();
		
		Arrays.fill(depths, -1);
		
		subV.clear();
		
		subV.reserve(G.numVertices());
		for (int v = 0; v < G.numVertices(); ++v)
		{
			if ((exclude != null) && exclude.isColored(v))
				continue;
			
			subV.pushBack(v);
		}
		assert !subV.empty();

		int n = subV.size();
		int skipped = subVC.prune(subV, lb, true);

		stack.add(new StackEntry(initBuffer, 0, skipped));
		
		if (verbose)
			logger.info(
				"Clique: t {}, 0+{}, V {}/{}/{}",
				t.now(), subVC.numColors(),
				subV.size()-skipped, skipped, n-subV.size());
		else
			logger.debug(
				"Clique: t {}, 0+{}, V {}/{}/{}",
				t.now(), subVC.numColors(),
				subV.size()-skipped, skipped, n-subV.size());
	}
	
	private int next(int lb, boolean verbose)
	{
		assert clique.size()+1 == stack.size();
		
		StackEntry entry = stack.getLast();

		// tried all
		if (entry.empty())
			return -1;
		
		if (entry.needUpdate())
		{
			int n = subV.size();
			int skipped = subVC.prune(subV,
				best.size()-clique.size()+1, true);
			entry.update(skipped);
			
			if (verbose)
				logger.info(
					"Clique Update: t {}, {}+{}, V {}/{}/{}",
					t.now(), clique.size(), subVC.numColors(),
					subV.size()-skipped, skipped, n-subV.size());
			else			
				logger.debug(
					"Clique Update: t {}, {}+{}, V {}/{}/{}",
					t.now(), clique.size(), subVC.numColors(),
					subV.size()-skipped, skipped, n-subV.size());

			// all pruned
			if (entry.empty())
				return -1;
		}
		
		return entry.popBack();
	}
	
	private void push(int v, int lb, boolean verbose)
	{
		assert clique.size()+1 == stack.size();

		subV.clear();
		for (int iv = G.neighborsBegin(v);
			iv != G.neighborsEnd(v); ++iv)
		{
			int u = G.neighborVertex(iv);
			
			// not at the same depth, ignore
			if (depths[u] != clique.size())
			{
				assert depths[u] < clique.size();
				continue;
			}
			
			subV.pushBack(u);
		}
		
		// pre-coloring check
		// prune less than lb or no more than best
		int max = subV.size()+1+clique.size();
		if ((max < lb) || (max <= best.size()))
			return;
		
		clique.pushBack(v);
		
		// update best
		if (best.size() < clique.size())
		{
			best.assign(clique);
			if (verbose)
				logger.info(
					"Clique Best: t {}, best {}, lb {}",
					t.now(), best.size(), lb);
			else
				logger.debug(
					"Clique Best: t {}, best {}, lb {}",
					t.now(), best.size(), lb);
		}
		
		if (subV.empty())
		{
			// no need to go further
			clique.popBack();
			return;
		}
		
		int n = subV.size();
		
		int lbEff = lb;
		if (lbEff < best.size()+1)
			lbEff = best.size()+1;
		
		int skipped = subVC.prune(subV,
			lbEff-clique.size(), false);
		
		if (subV.empty())
		{
			// no need to go further
			clique.popBack();
			return;
		}
		
		int lastSkipped = stack.getLast().skipped;
		int lastSize = stack.getLast().size;
		
		// one level deeper
		stack.add(new StackEntry(
			buffersArray, buffersOffsets[v], skipped));
		
		if (clique.size() == 1)
		{
			++l1Count;

			if ((l1Count%100 == 0) && verbose)
				logger.info(
					"Clique: t {}, {}+{}, V {}/{}/{}, {}/{}, push {}",
					t.now(), clique.size(), subVC.numColors(),
					subV.size()-skipped, skipped, n-subV.size(),
					lastSize-lastSkipped, lastSkipped, v);
			else
				logger.debug(
					"Clique: t {}, {}+{}, V {}/{}/{}, {}/{}, push {}",
					t.now(), clique.size(), subVC.numColors(),
					subV.size()-skipped, skipped, n-subV.size(),
					lastSize-lastSkipped, lastSkipped, v);
		}
		
	}
	
	private void pop()
	{
		assert clique.size()+1 == stack.size();
		
		// clear remaining
		stack.removeLast().clear();
		if (stack.isEmpty())
			return; // end of recursion
		
		clique.popBack();
	}
	
	private final GraphUV G;
	private final IntVector clique;
	private final IntVector best;
	
	private final VCPruner4 subVC;
	private final IntVector subV = new IntVector();
	
	private final int[] initBuffer;
	
	private final int[] buffersArray;
	private final int[] buffersOffsets;
	
	private class StackEntry
	{
		StackEntry(int[] buffer, int offset, int skipped)
		{
			this.buffer = buffer;
			this.offset = offset;
			
			update(skipped);
		}
		
		boolean needUpdate()
		{
			if (best == MaxClique4.this.best.size())
				return false;
			
			assert best < MaxClique4.this.best.size();
			
			subV.clear();
			for (int i = 0; i < size; ++i)
			{
				int v = buffer[offset+i];
				subV.pushBack(v);
				assert depths[v] == clique.size();
				--depths[v];
			}
			return true;
		}
		
		void update(int skipped)
		{
			for (int i = 0; i < subV.size(); ++i)
			{
				int v = subV.at(i);
				buffer[offset+i] = v;
				++depths[v];
				assert depths[v] == clique.size();
			}
			this.skipped = skipped;
			this.size = subV.size();
			this.best = MaxClique4.this.best.size();
		}
		
		boolean empty()
		{
			return size == skipped;
		}
		
		int popBack()
		{
			assert !empty();
			--size;
			int v = buffer[offset+size];
			assert depths[v] == clique.size();
			--depths[v];
			return v;
		}
		
		void clear()
		{
			for (int i = 0; i < size; ++i)
			{
				int v = buffer[offset+i];
				assert depths[v] == clique.size();
				--depths[v];
			}
			skipped = size = 0;
		}
		
		private final int offset;
		private final int[] buffer;

		private int skipped, size, best;
	}
	
	private final LinkedList<StackEntry> stack = new LinkedList<>();
	
	private final int[] depths;
	
	private final Timer t = new Timer();
	
	private int l1Count;
	
	private static final Logger logger
		= LoggerFactory.getLogger(MaxClique4.class);

}

class VCPruner4
{
	public VCPruner4(GraphUV G)
	{
		this.G = G;
		
		final int numV = G.numVertices();
		final int numE = G.numEdges();
		
		this.begins = new int[numV];
		this.ends = new int[numV];
		this.neighbors = new int[numE*2];
		
		this.colors = new int[numV];
		this.nbColors = new Coloring(numV);
		
		buffer.reserve(numV);
	}
	
	/**
	 * Color activeV to prune/order them.
	 * 
	 * @param activeV subgraph (in/out)
	 * @param lb prune any vertex not in k-clique 
	 * @param verbose show details
	 * @return # vertices to be skipped
	 */
	public int prune(IntVector activeV, int lb, boolean verbose)
	{
		final int n = activeV.size();
		if (lb < 1)
			lb = 1;
		
		// populate data structures
		init(activeV);
		
		int rounds = 1;
		for (;; ++rounds)
		{
			initColors(activeV);
			
			// generate coloring
			for (int i = 0; i < activeV.size(); ++i)
				colorOne(activeV.at(i));
			
			// prune vertices
			buffer.clear();
			for (int i = 0; i < activeV.size(); ++i)
			{
				int v = activeV.at(i);
				int m = countNeighborColors(v);
				
				// at most (m+1)-clique
				if (m+1 < lb)
				{
					// prune it and remove color
					// allow Gauss–Seidel iterations
					colors[v] = -1;
					continue;
				}
				
				buffer.pushBack(v);
			}
			
			if (buffer.size() == activeV.size())
				break;
			
			buffer.swap(activeV);
		}
		
		if (activeV.empty())
		{
			assert colorCounts.empty();
			
			if (verbose)
				logger.debug(
					"VCPrune: lb {}, rounds {}, {}->0",
					lb, rounds, n);
			return 0;
		}
		
		assert lb <= colorCounts.size();
		
		sortColors();
		for (int i = 0; i < activeV.size(); ++i)
			generateOne(activeV.at(i));
		buffer.swap(activeV);
		
		// skip vertices with max lb-1 colors
		int ret = offsets.at(lb-1);

		if (verbose)
			logger.debug(
				"VCPrune: lb {}, rounds {}, {}->{}/{}, {} {}",
				lb, rounds, n, activeV.size()-ret, ret,
				colorCounts.size(), colorCounts);
		
		return ret;
	}
	
	public int numColors()
	{
		return colorCounts.size();
	}

	private void init(IntVector activeV)
	{
		Coloring active = new Coloring(G.numVertices());		
		for (int i = 0; i < activeV.size(); ++i)
			active.setColor(activeV.at(i));
		
		for (int i = 0, e = 0; i < activeV.size(); ++i)
		{
			int v = activeV.at(i);
			
			begins[v] = e;
			
			for (int iv = G.neighborsBegin(v);
				iv != G.neighborsEnd(v); ++iv)
			{
				int u = G.neighborVertex(iv);
				if (active.isColored(u))
					neighbors[e++] = u;
			}
			
			ends[v] = e;
		}
		
		Integer[] sorted = new Integer[activeV.size()];
		for (int i = 0; i < activeV.size(); ++i)
			sorted[i] = activeV.at(i);
		
		Arrays.sort(sorted,
			(l, r) -> ends[r]-begins[r]-ends[l]+begins[l]);
		
		for (int i = 0; i < activeV.size(); ++i)
			activeV.set(i, sorted[i]);
	}
	
	private void initColors(IntVector activeV)
	{
		for (int i = 0; i < activeV.size(); ++i)
			colors[activeV.at(i)] = -1;
		
		colorCounts.clear();
	}
	
	private void colorOne(int v)
	{
		nbColors.resetAll();
		
		for (int iv = begins[v]; iv != ends[v]; ++iv)
		{
			int u = neighbors[iv];
			
			// no color, ignore
			if (colors[u] == -1)
				continue;
			
			nbColors.setColor(colors[u]);
		}
		
		int c = 0;
		for (; c < nbColors.size(); ++c)
		{
			if (!nbColors.isColored(c))
				break;
		}
		colors[v] = c;
		
		if (c == colorCounts.size())
		{
			colorCounts.pushBack(1);
		}
		else
		{
			assert c < colorCounts.size();
			int oldCount = colorCounts.at(c);
			colorCounts.set(c, oldCount+1);
		}
	}

	private int countNeighborColors(int v)
	{
		nbColors.resetAll();
		
		int ret = 0;
		for (int iv = begins[v]; iv != ends[v]; ++iv)
		{
			int u = neighbors[iv];
			
			// no color, ignore
			if (colors[u] == -1)
				continue;
			
			if (nbColors.isColored(colors[u]))
				continue;
			
			++ret;
			nbColors.setColor(colors[u]);
		}
		
		return ret;
	}
	
	private void sortColors()
	{
		Integer[] sorted = new Integer[colorCounts.size()];
		for (int c = 0; c < colorCounts.size(); ++c)
			sorted[c] = c;
		
		Arrays.sort(sorted,
			(l, r) -> colorCounts.at(r)-colorCounts.at(l));
		
		offsets.clear();
		toNewColor.assign(colorCounts.size(), -1);
		for (int nc = 0; nc < colorCounts.size(); ++nc)
		{
			int c = sorted[nc];
			
			toNewColor.set(c, nc);
			
			int prev = (nc == 0)? 0: offsets.back();
			
			offsets.pushBack(prev+colorCounts.at(c));
		}
		assert offsets.back() == buffer.size();
	}
	
	private void generateOne(int v)
	{
		int c = colors[v];
		int nc = toNewColor.at(c);
		
		int offset = offsets.at(nc)-1;
		
		buffer.set(offset, v);
		offsets.set(nc, offset);
	}
	
	private final GraphUV G;

	// vertices in subgraph
	private final int[] begins;
	private final int[] ends;
	private final int[] neighbors;
	
	// color per vertex
	private final int[] colors;

	// colors used by neighbors
	private final Coloring nbColors;
	
	// #vertices per color
	private final IntVector colorCounts = new IntVector();
	
	// mapped new color
	private final IntVector toNewColor = new IntVector();
	
	// offsets for counting sort the vertices
	private final IntVector offsets = new IntVector();
	
	// internal buffer
	private final IntVector buffer = new IntVector();
	
	private static final Logger logger
		= LoggerFactory.getLogger(VCPruner4.class);
}