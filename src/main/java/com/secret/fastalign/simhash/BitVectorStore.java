package com.secret.fastalign.simhash;

import com.secret.fastalign.data.Sequence;

public class BitVectorStore extends AbstractVectorStore<AbstractSequenceBitHash>
{

	public BitVectorStore(int kmerSize, int numWords) 
	{
		super(kmerSize, numWords);
	}

	@Override
	public AbstractSequenceBitHash getVectorHash(Sequence seq)
	{
		return new SequenceSimHash(seq, this.kmerSize, this.numWords);
	}

}
