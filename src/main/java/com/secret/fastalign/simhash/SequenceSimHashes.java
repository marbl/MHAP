package com.secret.fastalign.simhash;

import java.util.ArrayList;

import com.secret.fastalign.general.AbstractOrderedSequenceHashes;
import com.secret.fastalign.general.Sequence;

public final class SequenceSimHashes extends AbstractOrderedSequenceHashes<SimHash,SequenceSimHashes>
{
	public SequenceSimHashes(Sequence seq, int kmerSize, int numWords, int subStringSize, int subKmerSize, int subWordSize)
	{
		super(new SimHash(seq, kmerSize, numWords), generateSubHashes(seq, subStringSize, subKmerSize, subWordSize));
	}
	
	private static ArrayList<SimHash> generateSubHashes(Sequence seq, int subStringSize, int subKmerSize, int subWordSize)
	{
		//generate the array of simhashes
		ArrayList<SimHash> subHashes = new ArrayList<SimHash>(seq.length()-subKmerSize+1);
		for (int iter=0; iter<seq.length()-subKmerSize+1; iter+=subStringSize)
		{
			String subString = seq.getString().substring(iter, Math.min(iter+subStringSize, seq.length()));
			Sequence subSequence = new Sequence(subString, seq.getId());
			
			subHashes.add(new SimHash(subSequence, subKmerSize, subWordSize));
		}
		
		return subHashes;
	}
}
