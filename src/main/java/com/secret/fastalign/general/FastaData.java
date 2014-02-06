package com.secret.fastalign.general;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Utils;

public class FastaData
{
	// length of sequences loaded
	private final ConcurrentLinkedQueue<Sequence> sequenceList;
	
	private static final String[] fastaSuffix = {"fna", "contigs", "final", "fasta", "fa"};

	private FastaData(ConcurrentLinkedQueue<Sequence> seqList)
	{
		this.sequenceList = new ConcurrentLinkedQueue<Sequence>(seqList);
	}
	
	public FastaData(String file, int kmerSize) throws IOException 
	{
		BufferedReader bf;
		try
		{
			bf = Utils.getFile(file, fastaSuffix);
		}
		catch (Exception e)
		{
			throw new FastAlignRuntimeException(e);
		}
		
		String line = null;
		StringBuilder fastaSeq = new StringBuilder();

		this.sequenceList = new ConcurrentLinkedQueue<Sequence>();

		String header = "";
		while ((line = bf.readLine()) != null)
		{
			if (line.startsWith(">"))
			{
				if (fastaSeq.length() > 0)
					addMers(new SequenceId(header), fastaSeq.toString().toUpperCase());
				
				//header = line.substring(1).split("[\\s]+", 2)[0];
				
				//reset the storage
				fastaSeq.setLength(0);

			}
			else
			{
				fastaSeq.append(line);
			}
		}
		if (fastaSeq.length() != 0)
		{
			addMers(new SequenceId(header), fastaSeq.toString().toUpperCase());
		}
		bf.close();
	}
	
	// process a sequence and store the kmers/sequence length/etc
	public void addMers(SequenceId id, String seq)
	{
		Sequence sequence = new Sequence(seq, id);
		this.sequenceList.add(sequence);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public FastaData clone()
	{
		return new FastaData(this.sequenceList);
	}
	
	public Sequence dequeue()
	{
		return this.sequenceList.poll();
	}
	
	public Sequence getSequence(SequenceId id)
	{
		if (id.isForward())
		{
			for (Sequence seq : this.sequenceList)
				if (seq.getId().equals(id))
					return seq;
		}
	
		id = id.complimentId();
		for (Sequence seq : this.sequenceList)
			if (seq.getId().equals(id))
				return seq.getReverseCompliment();
	
		return null;
	}

	public ConcurrentLinkedQueue<Sequence> getSequences()
	{
		return this.sequenceList;
	}

	public boolean isEmpty()
	{
		return this.sequenceList.isEmpty();
	}

	public int size()
	{
		return this.sequenceList.size();
	}
}
