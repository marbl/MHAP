package com.secret.fastalign.simhash;

import com.secret.fastalign.data.Sequence;

public class KmerVectoreStore extends AbstractVectorStore<VectorSubKmerHash>
{

	public KmerVectoreStore(int kmerSize, int numWords)
	{
		super(kmerSize, numWords);
	}

	@Override
	public VectorSubKmerHash getVectorHash(Sequence seq)
	{
		return new VectorSubKmerHash(seq, this.kmerSize);
	}

}
