package com.secret.fastalign.hash;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;

import com.secret.fastalign.utils.Utils;

public class FastaData
{
	// length of sequences loaded
	private final ArrayList<Sequence> sequenceList;

	public FastaData(String file, String[] fastaSuffix, int kmerSize) throws Exception
	{
		BufferedReader bf = Utils.getFile(file, fastaSuffix);
		String line = null;
		StringBuilder fastaSeq = new StringBuilder();

		this.sequenceList = new ArrayList<Sequence>();

		// String header = "";
		long numProcessed = 0;
		while ((line = bf.readLine()) != null)
		{
			if (line.startsWith(">"))
			{
				if (fastaSeq.length() > 0)
					addMers(new SequenceId(numProcessed++), fastaSeq.toString().toUpperCase(), kmerSize);
				
				fastaSeq.setLength(0);
			}
			else
			{
				fastaSeq.append(line);
			}
		}
		if (fastaSeq.length() != 0)
		{
			addMers(new SequenceId(numProcessed++), fastaSeq.toString().toUpperCase(), kmerSize);
		}
		bf.close();

	}

	// process a sequence and store the kmers/sequence length/etc
	public void addMers(SequenceId id, String seq, int merSize)
	{
		this.sequenceList.add(new Sequence(seq, merSize, id));
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
