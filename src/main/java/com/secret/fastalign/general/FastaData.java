package com.secret.fastalign.general;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.secret.fastalign.utils.Utils;

public class FastaData
{
	// length of sequences loaded
	private final ArrayList<Sequence> sequenceList;
	private final HashMap<SequenceId, Sequence> sequenceMap;

	public FastaData(String file, String[] fastaSuffix, int kmerSize) throws Exception
	{
		BufferedReader bf = Utils.getFile(file, fastaSuffix);
		String line = null;
		StringBuilder fastaSeq = new StringBuilder();

		this.sequenceList = new ArrayList<Sequence>();
		this.sequenceMap = new HashMap<SequenceId, Sequence>();

		String header = "";
		while ((line = bf.readLine()) != null)
		{
			if (line.startsWith(">"))
			{
				if (fastaSeq.length() > 0)
					addMers(new SequenceId(header, true), fastaSeq.toString().toUpperCase(), kmerSize);
				
				header = line.substring(1).split("[\\s]+", 2)[0];
				
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
			addMers(new SequenceId(header, true), fastaSeq.toString().toUpperCase(), kmerSize);
		}
		bf.close();

	}

	// process a sequence and store the kmers/sequence length/etc
	public void addMers(SequenceId id, String seq, int merSize)
	{
		Sequence sequence = new Sequence(seq, id);
		this.sequenceList.add(sequence);
		this.sequenceMap.put(id, sequence);
	}
	
	public Sequence getSequence(int index)
	{
		return this.sequenceList.get(index);
	}
	
	public Sequence getSequence(SequenceId id)
	{
		if (id.isForward())
			return this.sequenceMap.get(id);

		Sequence seq = this.sequenceMap.get(id.complimentId());
		return seq.getReverseCompliment();
	}

	public List<Sequence> getSequences()
	{
		return this.sequenceList;
	}

	public int size()
	{
		return this.sequenceList.size();
	}
}
