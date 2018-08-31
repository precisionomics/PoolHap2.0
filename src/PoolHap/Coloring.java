package PoolHap; 

import java.util.Arrays;

/**
 * Fast coloring without the need to clear the memory
 * everytime before use.
 * <p>
 * Each vertex is either not colored or has a color of
 * 0 to <code>numColors()-1</code>.
 */
public class Coloring {

	private final int[] colors;
	
	// [0, ub) are valid, [0, lb) are not colored
	private int lb, ub;

	/**
	 * Initialize n vertices, each without a color,
	 * and prepare to use <code>numColors</code> colors.
	 * @param n
	 * @param numColors number of colors
	 */
	public Coloring(int n, int numColors)
	{
		this.colors = new int[n];
		this.lb = 1;
		this.ub = numColors+1;
	}
	
	/**
	 * Return number of vertices.
	 * @return number of vertices
	 */
	public int size() {return colors.length;}
	
	/**
	 * Return number of colors.
	 * @return number of colors
	 */
	public int numColors() {return ub-lb;}
	
	/**
	 * Reset all vertices to not colored
	 * and prepare to use <code>numColors</code> colors.
	 * @param numColors number of colors
	 */
	public void resetAll(int numColors)
	{
		assert numColors > 0;
		assert numColors+1 > 0;
		if (ub+numColors < 0) // overflow
		{
			lb = 1;
			ub = numColors+1;
			Arrays.fill(colors, 0);
		}
		else
		{
			lb = ub;
			ub = numColors+ub;
		}
	}
	
	/**
	 * Set the color of vertex <code>v</code> to <code>c</code>.
	 * @param v an vertex
	 * @param c the color
	 */
	public void setColor(int v, int c)
	{
		assert c >= 0;
		assert c < numColors();
		colors[v] = lb+c;
	}
	
	/**
	 * Reset the color of vertex <code>v</code> to be not colored.
	 * @param v an vertex
	 */
	public void resetColor(int v)
	{
		colors[v] = lb-1;
	}

	/**
	 * Check whether vertex <code>v</code> is colored or not
	 * @param v an vertex
	 * @return whether vertex <code>v</code> is colored or not
	 */
	public boolean isColored(int v)
	{
		return colors[v] >= lb;
	}
	
	/**
	 * Return the color of vertex <code>v</code> assuming it is colored.
	 * @param v an vertex with a color
	 * @return the color of vertex <code>v</code>
	 */
	public int getColor(int v)
	{
		assert isColored(v);
		return colors[v]-lb;
	}
	
	/**
	 * Initialize n vertices, each without a color,
	 * and prepare to use a single color.
	 * @param n
	 * <p>
	 * This is shortcut to {@link #Coloring(int, int) Coloring(n, 1)}. 
	 */
	public Coloring(int n)
	{
		this(n, 1);
	}
	
	/**
	 * Reset all vertices to not colored
	 * and prepare to use a single color.
	 * <p>
	 * This is shortcut to {@link #resetAll(int) resetAll(1)}. 
	 */
	public void resetAll() {resetAll(1);}
	
	/**
	 * Set vertex <code>v</code> to colored.
	 * <p>
	 * This is shortcut to {@link #setColor(int, int) setColor(v,0)}. 
	 * @param v an vertex
	 */
	public void setColor(int v) {setColor(v, 0);}

}
