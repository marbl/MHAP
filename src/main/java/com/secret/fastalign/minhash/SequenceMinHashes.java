package com.secret.fastalign.minhash;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.HashSet;

import com.secret.fastalign.general.OrderKmerHashes;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.general.SequenceHashes;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public final class SequenceMinHashes implements SequenceHashes
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3155689614837922443L;
	private final SequenceId id;
	private final MinHash mainHashes;
	private final OrderKmerHashes orderedHashes;
	
	public final static double SHIFT_CONSENSUS_PERCENTAGE = 0.75;
	
	public static SequenceMinHashes fromByteStream(DataInputStream input, int offset) throws IOException
	{
		try
		{
			input.
			
			//dos.writeInt(this.id.getHeaderId());
			//dos.writeBoolean(this.id.isForward());
			SequenceId id = new SequenceId(input.readInt()+offset, input.readBoolean());
			
			//dos.write(this.mainHashes.getAsByteArray());
			MinHash mainHashes = MinHash.fromByteStream(input);
			
			if (mainHashes==null)
				throw new FastAlignRuntimeException("Unexpected hash read error.");
			
			//dos.writeInt(this.completeHash.length);
			int hashLength = input.readInt();			
			
			int[][] completeHash = new int[hashLength][]; 
			for (int iter=0; iter<hashLength; iter++)
			{
				//dos.writeInt(this.completeHash[iter][iter2]);
				completeHash[iter] = new int[2];
				completeHash[iter][0] = input.readInt();
				completeHash[iter][1] = input.readInt();					
			}
			
			return new SequenceMinHashes(id, mainHashes, null);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	public SequenceMinHashes(SequenceId id, MinHash mainHashes, OrderKmerHashes orderedHashes)
	{
		this.id = id;
		this.mainHashes = mainHashes;
		this.orderedHashes = orderedHashes;
	}
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numHashes, int subSequenceSize, int orderedKmerSize, 
			boolean storeHashes, HashSet<Integer> filter)
	{
		this.id = seq.getId();
		this.mainHashes = new MinHash(seq, kmerSize, numHashes, subSequenceSize, filter);
		this.orderedHashes = new OrderKmerHashes(seq, orderedKmerSize);
	}
	
	public SequenceMinHashes createOffset(int offset)
	{
		return new SequenceMinHashes(this.id.createOffset(offset), this.mainHashes, this.orderedHashes);
	}
	
	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
    DataOutputStream dos = new DataOutputStream(bos);
    
    try
		{
      dos.writeInt(this.id.getHeaderId());
      dos.writeBoolean(this.id.isForward());
			dos.write(this.mainHashes.getAsByteArray());
			  
	    dos.flush();
	    return bos.toByteArray();
		}
    catch (IOException e)
    {
    	throw new FastAlignRuntimeException("Unexpected IO error.");
    }
	}

	public MinHash getMinHashes()
	{
		return this.mainHashes;
	}
	
	public OrderKmerHashes getOrderedHashes()
	{
		return this.orderedHashes;
	}
	
	@Override
	public SequenceId getSequenceId()
	{
		return this.id;
	}

	@Override
	public int getSequenceLength()
	{
		return this.mainHashes.getSequenceLength();
	}
}
