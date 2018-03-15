package PoolHap; 

import java.util.Arrays;
import java.util.Collection;

/*
 * There is no java equivalent of std::vector<int>.
 * Fortunately, it is no hard to implement one.
 */
public class IntVector {

	public IntVector() {}
	public IntVector(int n) {assign(n, 0);}
	public IntVector(int n, int val) {assign(n, val);}
	public IntVector(Collection<Integer> c) {assign(c);}
	public IntVector(int[] a) {assign(a);}
	public IntVector(int[] a, int b, int e) {assign(a, b, e);}
	public IntVector(IntVector v) {assign(v);}

	// capacity and size

	public void reserve(int n)
	{
		if (vec.length < n)
			vec = Arrays.copyOf(vec, n);
	}

	public void resize(int n, int val)
	{
		reserve(n);
		for (; last != n; ++last)
			vec[last] = val;
	}
	
	public int capacity() {return vec.length;}
	public void resize(int n) {reserve(last = n);}
	public void clear() {last = 0;}
	public int size() {return last;}
	public boolean empty() {return last == 0;}
	
	public void swap(IntVector rhs)
	{
		int l = last; last = rhs.last; rhs.last = l;
		int[] v = vec; vec = rhs.vec; rhs.vec = v;
	}

	// size of a vector never shrink unless told to do so
	
	public void assign(int n, int val)
	{
		reserve(n);
		Arrays.fill(vec, 0, last = n, val);
	}

	public void assign(Collection<Integer> c)
	{
		reserve(c.size());
		last = 0;
		for (Integer I: c)
			vec[last++] = I.intValue();
	}
	
	public void assign(int[] a)
	{
		reserve(a.length);
		System.arraycopy(a, 0, vec, 0, last = a.length);
	}

	public void assign(int[] a, int begin, int end)
	{
		reserve(end-begin);
		System.arraycopy(a, begin, vec, 0, last = end-begin);
	}

	public void assign(IntVector v)
	{
		reserve(v.size());
		System.arraycopy(v.data(), 0, vec, 0, last = v.size());
	}
	
	public int[] shrinkToFit()
	{
		if (last == vec.length)
			return vec;
		return vec = Arrays.copyOf(vec, last);
	}
	
	// accessors
	
	public int at(int i) {assert i < last; return vec[i];}
	public int get(int i) {return at(i);}
	public void set(int i, int val) {assert i < last; vec[i] = val;}
	
	public int front() {return at(0);}
	public int back() {return at(last-1);}
	public void setFront(int val) {set(0, val);}
	public void setBack(int val) {set(last-1, val);}
	
	// use with care, call shrink_to_fit() if necessary
	public int[] data() {return vec;}
	
	// modifiers
	
	public void insert(int i, int val)
	{
		insert(i, 1, val);
	}
	
	public void insert(int i, int n, int val)
	{
		assert i >= 0 && i <= last;
		assert n > 0;
		if (last+n < capacity())
			reserve(last+Math.max(n, last/2));
		last += n;
		for (int j = last-1; j != i+n-1; --j)
			vec[j] = vec[j-n];
		for (int k = 0; k < n; ++k)
			vec[i+k] = val;
	}
	
	public void erase(int i)
	{
		assert i >= 0 && i < last;
		for (--last; i != last; ++i)
			vec[i]= vec[i+1];
	}
	
	public int popBack()
	{
		assert last > 0;
		--last;
		return vec[last]; 
	}
	
	public void pushBack(int val)
	{
		if (last == capacity())
			reserve(last+Math.max(1, last/2));
		vec[last++] = val;
	}
	
	public void append(int a[])
	{
		if (last+a.length > capacity())
			reserve(last+Math.max(a.length, last/2));
		for (int i = 0; i < a.length; ++i)
			vec[last++] = a[i];
	}

	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append('[');
		if (size() != 0)
			sb.append(String.format("%d", vec[0]));
		for (int i = 1; i < size(); ++i)
		{
			sb.append(String.format(",%d", vec[i]));
		}
		sb.append(']');
		return sb.toString();
	}

	private int[] vec = new int[0];
	private int last = 0;

}
