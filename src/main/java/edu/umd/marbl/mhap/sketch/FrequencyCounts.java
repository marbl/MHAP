package edu.umd.marbl.mhap.sketch;

import java.util.HashMap;
import java.util.Map;


public final class FrequencyCounts
{
	private final double filterCutoff;
	private final Map<Long,Double> fractionCounts;
	private final double maxIdfValue;
	private final double maxValue;
	private final double minIdfValue;
	private final double minValue;
	
	public FrequencyCounts(Map<Long,Double> fractionCounts, double filterCutoff)
	{
		this.fractionCounts = new HashMap<>(fractionCounts);
		this.filterCutoff = filterCutoff;
		
		double maxValue = Double.NEGATIVE_INFINITY;
		for (double val : this.fractionCounts.values())
			maxValue = Math.max(maxValue, val);
		
		this.maxValue = maxValue;
		this.minValue = this.filterCutoff*0.1;
		
		this.minIdfValue = idf(this.maxValue);
		this.maxIdfValue = idf(this.minValue);
	}
	
	public boolean contains(long hash)
	{
		return this.fractionCounts.containsKey(hash);
	}
	
	public double documentFrequencyRatio(long hash)
	{
		Double val = this.fractionCounts.get(hash);
		if (val == null)
			val = this.minValue;
		
		return val;
	}
	
	public double getFilterCutoff()
	{
		return this.filterCutoff;
	}
	
	public double idf(double freq)
	{
		return Math.log(this.maxValue/freq);
		//return Math.log1p(this.maxValue/freq);
	}
	
	public double idf(long hash)
	{
		double freq = documentFrequencyRatio(hash);
		return idf(freq); 
	}
	
	public double idfDiscrete(long hash, int maxValue)
	{
		Double val = this.fractionCounts.get(hash);
		if (val == null)
			return maxValue;
		
		//get the true value
		double idf = idf(val);
		
		//scale it to match max
		double scale = (maxIdf()-minIdf())/(double)(maxValue-1.0);
		
		return 1.0+(idf-minIdf())/scale;
	}
	
	public double inverseDocumentFrequency(long hash)
	{
		return 1.0/documentFrequencyRatio(hash);
	}
	
	public double maxIdf()
	{
		return this.maxIdfValue;
	}
	
	public double minIdf()
	{
		return this.minIdfValue;
	}
}
