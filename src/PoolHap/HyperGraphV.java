package PoolHap; 

import java.util.Arrays;

// import misc.IntVector; Not necessary here because in the same project.

/**
 * Vector based hypergraph.
 */
public class HyperGraphV {

	public void clear()
	{
		offsets.clear();
		elements.clear();
		offsets.pushBack(0);
	}
	
	public void reserve(int nEdges, int nElements)
	{
		offsets.reserve(nEdges+1);
		elements.reserve(nElements);
	}
	
	public void shrinkToFit()
	{
		offsets.shrinkToFit();
		elements.shrinkToFit();
	}
	
	public int numEdges()
	{
		return offsets.size()-1;
	}
	
	public int numElements()
	{
		return elements.size();
	}
	
	public int numElements(int iEdge)
	{
		return offsets.at(iEdge+1)-offsets.at(iEdge);
	}

	public int getElement(int iEdge, int iElement)
	{
		assert iElement < numElements(iEdge);
		return elements.at(offsets.at(iEdge)+iElement);
	}
	
	public int[] elements(int iEdge)
	{
		return Arrays.copyOfRange(elements.data(),
			offsets.at(iEdge), offsets.at(iEdge+1));
	}

	public void appendEdge()
	{
		offsets.pushBack(elements.size());
	}

	public void appendElement(int element)
	{
		elements.pushBack(element);
		offsets.setBack(offsets.back()+1);
	}
	
	private IntVector offsets = new IntVector(1, 0);
	private IntVector elements = new IntVector();

}
