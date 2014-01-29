package com.secret.fastalign.simhash;

import java.util.ArrayList;
import com.secret.fastalign.data.Sequence;

public class BitVectorStore extends AbstractVectorStore<AbstractSequenceBitHash>
{
	public BitVectorStore(int kmerSize, int numWords) 
	{
		super(kmerSize, numWords);
	}

	@Override
	public SequenceVectorHashes<AbstractSequenceBitHash> getVectorHash(Sequence seq, int kmerSize, int numWords)
	{
		SequenceSimHash mainHash = new SequenceSimHash(seq, kmerSize, numWords);
		
		//generate the array of simhashes
		int len = seq.getString().length();
		ArrayList<SequenceSimHash> orderedHashes = new ArrayList<SequenceSimHash>(len-SUB_KMER_SIZE+1);
		for (int iter=0; iter<len-SUB_KMER_SIZE+1; iter+=SUB_STRING_SIZE)
		{
			String subString = seq.getString().substring(iter, Math.min(iter+SUB_STRING_SIZE, len));
			Sequence subSequence = new Sequence(subString, seq.getId());
			
			orderedHashes.add(new SequenceSimHash(subSequence, SUB_KMER_SIZE, SUB_WORD_SIZE));
		}
		
		return new SequenceVectorHashes<AbstractSequenceBitHash>(len, mainHash, orderedHashes);
	}

}
