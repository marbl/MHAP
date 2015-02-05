package edu.umd.marbl.mhap.align;

import java.util.Iterator;
import java.util.List;

public final class Alignment<S extends AlignElement<S>>
{
	public enum Operation
	{
		MATCH,
		INSERT,
		DELETE;
	}
	
	private final double score;
	//private final double gapOpen;
	private final List<Operation> operations;
	private final S a;
	private final S b;
	
	protected Alignment(S a, S b, double score, double gapOpen, List<Operation> operations)
	{
		this.score = score;
		this.operations = operations;
		this.a = a;
		this.b = b;
		//this.gapOpen = gapOpen;
	}
	
	public double getOverlapScore(int minMatches)
	{
		int i = 0;
		int j = 0;
		
		Iterator<Operation> iter = this.operations.iterator();
		
		if (!iter.hasNext())
			return 0.0;
		
		//remove the start
		Operation o = iter.next();			
		while (o==Operation.DELETE)
		{
			i++;
			
			if (iter.hasNext())
				o = iter.next();
			else
				return 0.0;
		}
		if (i==0)
		{
			while (o==Operation.INSERT)
			{
				if (iter.hasNext())
					o = iter.next();
				else
					return 0.0;
			}
		}
		
		double score = 0.0;
		int count = 0;
		while (true)
		{
			if (o == Operation.DELETE)
			{	
				i++;
			}
			else
			if (o == Operation.INSERT)
			{
				//count++;
				j++;
			}
			else
			{
				score = score + a.similarityScore(b, i, j)-a.getSimOffset();
				count++;
				i++;
				j++;
			}
			
			if (iter.hasNext())
				o = iter.next();
			else
				break;
		}
		
		System.err.println(this.operations);
		System.err.println("HI="+count+" "+minMatches);
		
		if (count<minMatches)
			return 0.0;
		
		if (score<=0.0)
			return 0.0;
		
		return score/(double)(count);
	}
	
	public double getScore()
	{
		return this.score;
	}
	
	public int getA1()
	{
		int start = 0;
				
		return start;
	}
	
	public int getA2()
	{
		int start = 0;
				
		return start;
	}
	
	public int getB1()
	{
		int end = 0;
				
		return this.a.length()-end;
	}
	
	public int getB2()
	{
		int end = 0;
				
		return this.b.length()-end;
	}
	
	public String outputAlignmentSelf()
	{
		StringBuilder str = new StringBuilder();

		int i = 0;
		int j = 0;
		
		int count = 0;
		while(i<a.length() || j<b.length())
		{
			Operation o;
			if (count<this.operations.size())
				o = this.operations.get(count);
			else
			if (i<a.length())
				o = Operation.DELETE;
			else
				o = Operation.INSERT;
						
			if (o == Operation.DELETE)
			{	
				String aStr;
				if (j>=b.length())
					aStr = a.toString(i);
				else
					aStr = a.toString(b, i, j);
				
				str.append(aStr);
				i++;
			}
			else
			if (o == Operation.INSERT)
			{
				String bStr;
				if (i>=a.length())
					bStr = b.toString(j);
				else
					bStr = b.toString(a, j, i);
				for (int space=0; space<bStr.length(); space++)
					str.append("-");
				j++;
			}
			else
			{
				str.append(a.toString(b, i, j));
				i++;
				j++;
			}
			
			count++;
		}
		
		return str.toString();		
	}
	
	public String outputAlignmentOther()
	{
		StringBuilder str = new StringBuilder();

		int i = 0;
		int j = 0;
		int count = 0;
		while(i<a.length() || j<b.length())
		{
			Operation o;
			if (count<this.operations.size())
				o = this.operations.get(count);
			else
			if (i<a.length())
				o = Operation.DELETE;
			else
				o = Operation.INSERT;
			
			if (o == Operation.INSERT)
			{
				String bStr;
				if (i>=a.length())
					bStr = b.toString(j);
				else
					bStr = b.toString(a, j, i);

				str.append(bStr);
				j++;
			}
			else
			if (o == Operation.DELETE)
			{
				String aStr;
				if (j>=b.length())
					aStr = a.toString(i);
				else
					aStr = a.toString(b, i, j);
				
				for (int space=0; space<aStr.length(); space++)
					str.append("-");
				i++;
			}
			else
			{
				str.append(b.toString(a, j, i));
				i++;
				j++;
			}
			
			count++;
		}
		
		return str.toString();
	}
	
	public String outputAlignment()
	{
		StringBuilder str = new StringBuilder();
		
		str.append("Sequence 1: ");
		str.append(outputAlignmentSelf()+"\n");
		str.append("Sequence 2: ");
		str.append(outputAlignmentOther()+"\n");
		
		return str.toString();
	}
}
