package com.secret.fastalign.direct;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import com.secret.fastalign.general.AbstractSequenceHashStreamer;
import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.ReadBuffer;
import com.secret.fastalign.utils.Utils;

public class SequenceDirectHashStreamer extends AbstractSequenceHashStreamer<SequenceDirectHashes>
{
	private final DataInputStream buffInput;
	private final HashSet<Integer> filter;
	private final int kmerSize;
	private final int offset;
	private boolean readClosed;
	private final int orderedKmerSize;

	public SequenceDirectHashStreamer(String file, int offset) throws FileNotFoundException
	{
		super(null, false);
		
		this.kmerSize = 0;
		this.orderedKmerSize = 0;
		this.filter = null;
		this.readClosed = false;
		this.offset = offset;

		this.buffInput = new DataInputStream(new BufferedInputStream(new FileInputStream(file), Utils.BUFFER_BYTE_SIZE));  	
	}
	
	public SequenceDirectHashStreamer(String file, int kmerSize, int orderedKmerSize,
			HashSet<Integer> filter, int offset) throws IOException
	{	
		super(new FastaData(file, offset), true);

		this.kmerSize = kmerSize;
		this.orderedKmerSize = orderedKmerSize;
		this.filter = filter;
		this.buffInput = null;
		this.readClosed = false;
		this.offset = offset;
	}

	@Override
	public SequenceDirectHashes getHashes(Sequence seq)
	{
		//compute the hashes
		return new SequenceDirectHashes(seq, this.kmerSize, this.orderedKmerSize, this.filter);
	}
	
	@Override
	protected SequenceDirectHashes readFromBinary(ReadBuffer buf, boolean fwdOnly) throws IOException
	{
		byte[] byteArray = null;
		synchronized (this.buffInput)
		{
			if (this.readClosed)
				return null;
			
			try
			{
				boolean keepReading = true;
				while(keepReading)
				{
					byte isFwd = this.buffInput.readByte();
					
					if (!fwdOnly || isFwd==1)
						keepReading = false;
					
					//get the size in bytes
					int byteSize = this.buffInput.readInt();
					
					//allocate the array
					byteArray = buf.getBuffer(byteSize);
					//byteArray = new byte[byteSize];				
					
					//read that many bytes
					this.buffInput.read(byteArray, 0, byteSize);
				}
			}
			catch(EOFException e)
			{
	  		this.buffInput.close();
	  		this.readClosed = true;			
	  		
	  		return null;
			}
		}
			
		//get as byte array stream
		SequenceDirectHashes seqHashes = SequenceDirectHashes.fromByteStream(new DataInputStream(new ByteArrayInputStream(byteArray)), this.offset);

  	return seqHashes;
  }
}
