package PoolHap; 

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 * Vector based read-only undirected graph.
 */

public class GraphUV {

	/**
	 * Neighbors are guaranteed to be stable as edges.
	 */
	public GraphUV(int numV, int[] edges)
	{
		this.numV = numV;
		this.edges = edges;
		this.nbBase = new int[numV+1];
		this.nbVertices = new int[edges.length];
		this.nbEdges = new int[edges.length];
		
		initializeNeighbors();
	}
	
	public int numVertices()
	{
		return numV;
	}
	
	public int numEdges()
	{
		return edges.length/2;
	}
	
	public int fromVertex(int e)
	{
		return edges[e*2];
	}
	public int toVertex(int e)
	{
		return edges[e*2+1];
	}
	
	public int numNeighbors(int v)
	{
		return nbBase[v+1]-nbBase[v];
	}
	
	public int neighborsBegin(int v)
	{
		return nbBase[v];
	}
	public int neighborsEnd(int v)
	{
		return nbBase[v+1];
	}
	
	public int neighborVertex(int i)
	{
		return nbVertices[i];
	}
	public int neighborEdge(int i)
	{
		return nbEdges[i];
	}
	
	public int firstNeighborVertex(int v)
	{
		assert numNeighbors(v) > 0;
		return neighborVertex(neighborsBegin(v));
	}
	public int lastNeighborVertex(int v)
	{
		assert numNeighbors(v) > 0;
		return neighborVertex(neighborsEnd(v)-1);
	}
	
	public int[] neighborVertices(int v)
	{
		return Arrays.copyOfRange(nbVertices,
			nbBase[v], nbBase[v+1]);
	}
	public int[] neighborEdges(int v)
	{
		return Arrays.copyOfRange(nbEdges,
			nbBase[v], nbBase[v+1]);
	}
	
	private void initializeNeighbors()
	{
		// count vertices' degree
		int[] degree = new int[numV];
		
		for (int e = 0; e < numEdges(); ++e)
		{
			++degree[toVertex(e)];
			++degree[fromVertex(e)];
		}

		// make adjacent index
		nbBase[0] = 0;
		for (int v = 1; v < numV; ++v)
			nbBase[v] = nbBase[v-1]+degree[v-1];
		nbBase[numV] = numEdges()*2;

		// fill adjacent list
		for (int e = 0; e < numEdges(); ++e)
		{
			int from = fromVertex(e);
			int to = toVertex(e);
			
			// fill from's neighbor: stable ordering
			int outi = nbBase[from+1]-degree[from];
			nbEdges[outi] = e;
			nbVertices[outi] = to;
			--degree[from];

			// fill to's neighbor: stable ordering
			int ini = nbBase[to+1]-degree[to];
			nbEdges[ini] = e;
			nbVertices[ini] = from;
			--degree[to];
		}
	}

	private final int numV;
	private final int[] edges;
	
	private final int[] nbBase; 
	private final int[] nbVertices;
	private final int[] nbEdges;
	
	/**
	 * Save the graph in the DIMACS format
	 * @param file the file to be saved
	 */
	public void saveDIMACS(File file)
		throws Exception
	{
		try(
			FileOutputStream f = new FileOutputStream(file);
			OutputStreamWriter osr = new OutputStreamWriter(f, "UTF-8");
			PrintWriter pw = new PrintWriter(osr))
		{
			pw.printf("p %s %d %d%n",
				file.getName(), numVertices(), numEdges());
			
			for (int e = 0; e < numEdges(); ++e)
				pw.printf("e %d %d%n",
					fromVertex(e)+1, toVertex(e)+1);
		}
	}	
	
	/**
	 * Apply BFS to discover connected components.
	 * @return the vertices in each component.
	 */
	public HyperGraphV bfsCC()
	{
		HyperGraphV ccs = new HyperGraphV();
		ccs.reserve(0, numVertices());
		
		Coloring colors = new Coloring(numVertices());
		
		int[] Q = new int[numVertices()];
		
		for (int u = 0; u < numVertices(); ++u)
		{
			if (colors.isColored(u))
				continue;
			
			// new cc
			ccs.appendEdge();
			
			// BFS from u
			Q[0] = u;
			colors.setColor(u);
			for (int begin = 0, end = 1;
				begin != end; ++begin)
			{
				// new vertex in the cc
				int v = Q[begin];
				ccs.appendElement(v);
				
				// check neighbors
				for (int i = neighborsBegin(v);
					i != neighborsEnd(v); ++i)
				{
					int nb = neighborVertex(i);
					if (colors.isColored(nb))
						continue; // saw it before
					
					// queue it
					Q[end++] = nb;
					colors.setColor(nb);
				}
			}
		}
		
		return ccs;
	}
	
	public static class SubGraph extends GraphUV
	{
		private SubGraph(GraphUV parent,
			int[] subEdges, int[] parentV, int[] parentE)
		{
			super(parentV.length, subEdges);
			
			this.parent = parent;
			this.parentV = parentV;
			this.parentE = parentE;
		}
		
		public GraphUV getParent()
		{
			return parent;
		}
		
		public int parentVertex(int v)
		{
			return parentV[v];
		}
		
		public int parentEdge(int e)
		{
			return parentE[e];
		}
		
		private final GraphUV parent;
		private final int[] parentV;
		private final int[] parentE;
	}
	
	/**
	 * Create a subgraph from a set of vertices.
	 * @param parentV vertices to be kept
	 * @return the subgraph
	 */
	public SubGraph project(int[] parentV)
	{
		Arrays.sort(parentV);
		
		// map from parent to sub
		int[] vMap = new int[numVertices()];
		Arrays.fill(vMap, -1);
		for (int subV = 0; subV < parentV.length; ++subV)
			vMap[parentV[subV]] = subV;
		
		IntVector subEdges = new IntVector();
		IntVector parentE = new IntVector();
		subEdges.reserve(numEdges()*2); // pairs
		parentE.reserve(numEdges());
		
		for (int e = 0; e < numEdges(); ++e)
		{
			int subFrom = vMap[fromVertex(e)];
			int subTo = vMap[toVertex(e)];
			
			if ((subFrom == -1) || (subTo == -1))
				continue; // not in the subgraph
			
			subEdges.pushBack(subFrom);
			subEdges.pushBack(subTo);
			parentE.pushBack(e);
		}
		
		return new SubGraph(this, subEdges.shrinkToFit(),
			parentV, parentE.shrinkToFit());
	}
	
	/**
	 * Decompose into a set of subgraphs using sets of vertices.
	 * @param sets sets of subgraph vertices
	 * @return sets of subgraphs
	 */
	public SubGraph[] decompose(HyperGraphV sets)
	{
		final int n = sets.numEdges();
		
		int[] subMap = new int[numVertices()];
		int[] vMap = new int[numVertices()];

		// -1 for vertices not belong to any subgraph.
		Arrays.fill(subMap, -1);
		
		// for stable ordering
		int[][] parentV = new int[n][];
		
		// edges
		IntVector[] subEdges = new IntVector[n];
		IntVector[] parentE = new IntVector[n];

		for (int sub = 0; sub < n; ++sub)
		{
			parentV[sub] = sets.elements(sub);
			
			// stable ordering
			Arrays.sort(parentV[sub]);
			
			for (int subV = 0;
				subV < parentV[sub].length;
				++subV)
			{
				int v = parentV[sub][subV];
				subMap[v] = sub;
				vMap[v] = subV;
			}

			subEdges[sub] = new IntVector();
			parentE[sub] = new IntVector();
		}
		
		for (int e = 0; e < numEdges(); ++e)
		{
			int from = fromVertex(e);
			int to = toVertex(e);
			
			if ((subMap[from] != subMap[to])
				|| (subMap[from] == -1))
				continue; // not a subgraph edge
			
			int sub = subMap[from];
			int subFrom = vMap[from];
			int subTo = vMap[to];
			
			subEdges[sub].pushBack(subFrom);
			subEdges[sub].pushBack(subTo);
			parentE[sub].pushBack(e);
		}
		
		SubGraph[] subs = new SubGraph[n];
		for (int sub = 0; sub < n; ++sub)
			subs[sub] = new SubGraph(this,
				subEdges[sub].shrinkToFit(),
				parentV[sub],
				parentE[sub].shrinkToFit());
		
		return subs;
	}

}
