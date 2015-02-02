package edu.umd.marbl.mhap.align;

import java.util.ArrayList;

public final class Alignment<S extends AlignElement<S>>
{
	public enum Operation
	{
		MATCH,
		INSERT,
		DELETE;
	}
	
	private final double score;
	private final double gapOpen;
	private final ArrayList<Operation> operations;
	private final S a;
	private final S b;
	
	protected Alignment(S a, S b, double score, double gapOpen, ArrayList<Operation> operations)
	{
		this.score = score;
		this.operations = operations;
		this.a = a;
		this.b = b;
		this.gapOpen = gapOpen;
		
		System.err.println(operations);
	}
	
	public double getOverlapScore()
	{
		int start = 0;
		int end = operations.size()-1;
		
		if (this.operations.get(start) == Operation.INSERT)
			while (start<=end && this.operations.get(start) == Operation.INSERT)
				start++;
		else
		if (this.operations.get(end) == Operation.DELETE)
			while (start<=end && this.operations.get(start) == Operation.DELETE)
				start++;
		
		if (this.operations.get(end) == Operation.INSERT)
			while (start<=end &&  this.operations.get(end) == Operation.INSERT)
				end--;
		else
		if (this.operations.get(end) == Operation.DELETE)
			while (start<=end &&  this.operations.get(end) == Operation.DELETE)
				end--;
		
		int i = 0;
		int j = 0;

		double score = 0.0;
		for (int iter=end; iter>=start; iter--)
		{
			Operation o = this.operations.get(iter);
			if (o == Operation.INSERT)
			{	
				score = score+gapOpen;
				i++;
			}
			else
			if (o == Operation.DELETE)
			{
				score = score+gapOpen;
				j++;
			}
			else
			{
				score = score + a.similarityScore(b, i, j);
				i++;
				j++;
			}
		}
		
		if (score<=0.0)
			return 0.0;
		
		return score/(double)(end-start+1);
	}
	
	public double getScore()
	{
		return this.score;
	}
	
	public int getA1()
	{
		int index = operations.size()-1;
		int start = 0;
		
		if (this.operations.get(index) == Operation.INSERT)
		{
			while (index>=0 && this.operations.get(index) == Operation.INSERT)
			{
				start++;
				index--;
			}
		}
		
		return start;
	}
	
	public int getA2()
	{
		int index = operations.size()-1;
		int start = 0;
		
		if (this.operations.get(index) == Operation.DELETE)
		{
			while (index>=0 && this.operations.get(index) == Operation.DELETE)
			{
				start++;
				index--;
			}
		}
		
		return start;
	}
	
	public int getB1()
	{
		int index = 0;
		int end = 0;
		
		if (this.operations.get(index) == Operation.INSERT)
		{
			while (index<this.operations.size() && this.operations.get(index) == Operation.INSERT)
			{
				end++;
				index++;
			}
		}
		
		return this.a.length()-end;
	}
	
	public int getB2()
	{
		int index = 0;
		int end = 0;
		
		if (this.operations.get(index) == Operation.DELETE)
		{
			while (index<this.operations.size() && this.operations.get(index) == Operation.DELETE)
			{
				end++;
				index++;
			}
		}
		
		return this.b.length()-end;
	}
	
	public String outputAlignmentSelf()
	{
		StringBuilder str = new StringBuilder();

		int i = 0;
		int j = 0;
		
		for (int iter=operations.size()-1; iter>=0; iter--)
		{
			Operation o = this.operations.get(iter);
			if (o == Operation.INSERT)
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
			if (o == Operation.DELETE)
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
		}
		
		return str.toString();		
	}
	
	public String outputAlignmentOther()
	{
		StringBuilder str = new StringBuilder();

		int i = 0;
		int j = 0;
		for (int iter=operations.size()-1; iter>=0; iter--)
		{
			Operation o = this.operations.get(iter);
			if (o == Operation.DELETE)
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
			if (o == Operation.INSERT)
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
