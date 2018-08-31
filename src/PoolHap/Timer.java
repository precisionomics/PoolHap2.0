package PoolHap; 

public class Timer {

	public Timer()
	{
		this.start = System.currentTimeMillis();
	}
	
	public double now()
	{
		return (System.currentTimeMillis()-start)/1000.0;
	}
	
	private final long start;

}
