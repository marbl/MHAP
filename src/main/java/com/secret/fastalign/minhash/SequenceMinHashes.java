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
			//input.
			
			//dos.writeBoolean(this.id.isForward());
			boolean isFwd = input.readBoolean();
			
			//dos.writeInt(this.id.getHeaderId());
			SequenceId id = new SequenceId(input.readInt()+offset, isFwd);
			
			//dos.write(this.mainHashes.getAsByteArray());
			MinHash mainHashes = MinHash.fromByteStream(input);
			
			if (mainHashes==null)
				throw new FastAlignRuntimeException("Unexpected data read error.");
			
			//dos.write(this.orderedHashes.getAsByteArray());
			OrderKmerHashes orderedHashes = OrderKmerHashes.fromByteStream(input);
			if (orderedHashes==null)
				throw new FastAlignRuntimeException("Unexpected data read error.");
			
			return new SequenceMinHashes(id, mainHashes, orderedHashes);
			
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
	
	@Override
	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream(this.mainHashes.numHashes()*4+this.orderedHashes.size()*2);
    DataOutputStream dos = new DataOutputStream(bos);
    
    try
		{
      dos.writeBoolean(this.id.isForward());
      dos.writeInt(this.id.getHeaderId());
			dos.write(this.mainHashes.getAsByteArray());
			dos.write(this.orderedHashes.getAsByteArray());
			
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
