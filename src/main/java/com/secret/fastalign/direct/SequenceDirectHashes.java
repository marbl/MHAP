package com.secret.fastalign.direct;

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

public class SequenceDirectHashes implements SequenceHashes
{

	/**
	 * 
	 */
	private static final long serialVersionUID = -1524284712894301880L;

	private final SequenceId id;
	private final OrderKmerHashes mainHashes;
	private final OrderKmerHashes orderedHashes;
	private final int seqLength;
	
	public final static double SHIFT_CONSENSUS_PERCENTAGE = 0.75;
	
	public static SequenceDirectHashes fromByteStream(DataInputStream input, int offset) throws IOException
	{
		try
		{
			//input.
			
			//dos.writeBoolean(this.id.isForward());
			boolean isFwd = input.readBoolean();
			
			//dos.writeInt(this.id.getHeaderId());
			SequenceId id = new SequenceId(input.readInt()+offset, isFwd);
			
			//dos.write(this.orderedHashes.getAsByteArray());
			OrderKmerHashes orderedHashes = OrderKmerHashes.fromByteStream(input);
			if (orderedHashes==null)
				throw new FastAlignRuntimeException("Unexpected data read error.");
			
			return new SequenceDirectHashes(id, 0, orderedHashes);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	public SequenceDirectHashes(SequenceId id, int seqLength, OrderKmerHashes orderedHashes)
	{
		this.id = id;
		this.seqLength = seqLength;
		this.orderedHashes = orderedHashes;
		this.mainHashes = null; //TODO fix
	}
	
	public SequenceDirectHashes(Sequence seq, int kmerSize, int orderedKmerSize, HashSet<Integer> filter)
	{
		this.id = seq.getId();
		this.seqLength = seq.length();
		this.mainHashes = new OrderKmerHashes(seq, kmerSize);
		this.orderedHashes = new OrderKmerHashes(seq, orderedKmerSize);
	}
	
	public SequenceDirectHashes createOffset(int offset)
	{
		return new SequenceDirectHashes(this.id.createOffset(offset), this.seqLength, this.orderedHashes);
	}
	
	@Override
	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream(this.orderedHashes.size()*2);
    DataOutputStream dos = new DataOutputStream(bos);
    
    try
		{
      dos.writeBoolean(this.id.isForward());
      dos.writeInt(this.id.getHeaderId());
			dos.write(this.orderedHashes.getAsByteArray());
			
	    dos.flush();
	    return bos.toByteArray();
		}
    catch (IOException e)
    {
    	throw new FastAlignRuntimeException("Unexpected IO error.");
    }
	}

	public OrderKmerHashes getMainHashes()
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
		return this.seqLength;
	}

}
