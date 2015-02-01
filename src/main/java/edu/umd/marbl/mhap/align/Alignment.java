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
	private final ArrayList<Operation> operations;
	private final S a;
	private final S b;
	
	protected Alignment(S a, S b, double score, ArrayList<Operation> operations)
	{
		this.score = score;
		this.operations = operations;
		this.a = a;
		this.b = b;
	}
	
	public double getScore()
	{
		return this.score;
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
