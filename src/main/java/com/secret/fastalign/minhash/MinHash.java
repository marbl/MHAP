package com.secret.fastalign.minhash;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashSet;

import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Utils;

public final class MinHash
{
	private final int[][] minHashes;
	private final int seqLength;
	
	public static MinHash fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			//store the size
			//bb.putInt(this.seqLength);
			int seqLength = input.readInt();
			
			//bb.putInt(this.minHashes.length);
			int seqNum = input.readInt();
			
			//bb.putInt(this.minHashes[0].length);
			int hashNum = input.readInt();
			
			//store the array
			int[][] minHashes = new int[seqNum][];
			for (int seq=0; seq<seqNum; seq++)
			{
				minHashes[seq] = new int[hashNum];
				for (int hash=0; hash<hashNum; hash++)
				{
					//bb.putInt(this.minHashes[seq][hash]);
					minHashes[seq][hash] = input.readInt();
				}
			}
			
			return new MinHash(seqLength, minHashes);
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	private MinHash(int seqLength, int[][] minHashes)
	{
		this.seqLength = seqLength;
		this.minHashes = minHashes;
	}
	
	public MinHash(Sequence seq, int kmerSize, int numHashes, int subSequenceSize, HashSet<Integer> filter)
	{
		this.seqLength = seq.length();

		int numberSubSeq = seq.length()/subSequenceSize+1;
		if (seq.length()%subSequenceSize<kmerSize)
			numberSubSeq--;
		
		//adjust for more equal distribution
		subSequenceSize = seq.length()/numberSubSeq+1;
		
		this.minHashes = new int[numberSubSeq][];
		for (int iter=0; iter<numberSubSeq; iter++)
		{
			String subString = seq.getString().substring(iter*subSequenceSize, Math.min(seq.length(), (iter+1)*subSequenceSize));

			//get the hashes
			if (subString.length()>=kmerSize)
				this.minHashes[iter] = Utils.computeKmerMinHashes(subString, kmerSize, numHashes, filter);
		}	
	}

	public byte[] getAsByteArray()
	{
		ByteBuffer bb = ByteBuffer.allocate(4*(3+this.minHashes.length*this.minHashes[0].length));
		
		//store the size
		bb.putInt(this.seqLength);
		bb.putInt(this.minHashes.length);
		bb.putInt(this.minHashes[0].length);
		
		//store the array
		for (int seq=0; seq<this.minHashes.length; seq++)
			for (int hash=0; hash<this.minHashes[0].length; hash++)
				bb.putInt(this.minHashes[seq][hash]); 
    
    return bb.array();
	}
	
	public final int getSequenceLength()
	{
		return this.seqLength;
	}

	/**
	 * @return the minHashes
	 */
	public final int[][] getSubSeqMinHashes()
	{
		return this.minHashes;
	}
	
	public final int numSubSequences()
	{
		return this.minHashes.length;
	}
	
	public final int numHashes()
	{
		return this.minHashes[0].length;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "MinHash "+Arrays.toString(this.minHashes) + "";
	}
}
