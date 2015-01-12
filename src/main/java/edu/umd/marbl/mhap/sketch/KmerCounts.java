package edu.umd.marbl.mhap.sketch;

public class KmerCounts
{
	private final CountMin<Long> counter;
	private final long totalReads;
	
	public KmerCounts(CountMin<Long> counter, int totalReads)
	{
		this.counter = counter;
		this.totalReads = totalReads;
	}
	
	public long getTotalReads()
	{
		return this.totalReads;
	}
	
	public double inverseDocumentFrequency(long kmer)
	{
		long count = counter.getCount(kmer);
		long total = counter.totalAdded();
		
		return Math.log10((double)total/(double)(count+1)); 
	}
	
	public double documentFrequencyRatio(long kmer)
	{
		long count = counter.getCount(kmer);
		long total = counter.totalAdded();
		
		return (double)count/(double)total; 
	}
	
	public double weight(long kmer, int countInRead, int numKmers)
	{
		double freqInData = inverseDocumentFrequency(kmer);
		
		//double freqInRead = (double)countInRead/(double)numKmers;
		
		//bias adjust
		double freqInRead = 0.5+0.5*(double)countInRead/numKmers;
		
		return freqInRead*freqInData; 
	}
}
