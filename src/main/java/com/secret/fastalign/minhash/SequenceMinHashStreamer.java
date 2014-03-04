package com.secret.fastalign.minhash;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.general.AbstractSequenceHashStreamer;
import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.ReadBuffer;
import com.secret.fastalign.utils.Utils;

public class SequenceMinHashStreamer extends AbstractSequenceHashStreamer<SequenceMinHashes>
{
	private final DataInputStream buffInput;
	private final HashSet<Integer> filter;
	private final int kmerSize;
	private final AtomicLong numberSubSequencesProcessed;
	private final int numHashes;
	private final int offset;
	private boolean readClosed;
	private final int orderedKmerSize;
	private final int subSequenceSize;

	public SequenceMinHashStreamer(String file, int offset) throws FileNotFoundException
	{
		super(null, false);
		
		this.kmerSize = 0;
		this.numHashes = 0;
		this.subSequenceSize = 0;
		this.orderedKmerSize = 0;
		this.filter = null;
		this.numberSubSequencesProcessed = new AtomicLong();
		this.readClosed = false;
		this.offset = offset;

		this.buffInput = new DataInputStream(new BufferedInputStream(new FileInputStream(file), Utils.BUFFER_BYTE_SIZE));  	
	}
	
	public SequenceMinHashStreamer(String file, int kmerSize, int numHashes, int subSequenceSize, int orderedKmerSize, 
			HashSet<Integer> filter, int offset) throws IOException
	{	
		super(new FastaData(file, offset), true);

		this.kmerSize = kmerSize;
		this.numHashes = numHashes;
		this.subSequenceSize = subSequenceSize;
		this.orderedKmerSize = orderedKmerSize;
		this.filter = filter;
		this.numberSubSequencesProcessed = new AtomicLong();
		this.buffInput = null;
		this.readClosed = false;
		this.offset = offset;
	}

	@Override
	public SequenceMinHashes getHashes(Sequence seq)
	{
		//compute the hashes
		return new SequenceMinHashes(seq, this.kmerSize, this.numHashes, this.subSequenceSize, this.orderedKmerSize, false, this.filter);
	}
	
	@Override
	public int getNumberSubSequencesProcessed()
	{
		return this.numberSubSequencesProcessed.intValue();
	}

	@Override
	protected void processAddition(SequenceMinHashes seqHashes)
	{
		super.processAddition(seqHashes);
		if (seqHashes!=null)
			this.numberSubSequencesProcessed.getAndAdd(seqHashes.getMinHashes().numSubSequences());
	}

	@Override
	protected SequenceMinHashes readFromBinary(ReadBuffer buf) throws IOException
	{
		byte[] byteArray = null;
		synchronized (this.buffInput)
		{
			if (this.readClosed)
				return null;
			
			try
			{
				//get the size in bytes
				int byteSize = this.buffInput.readInt();
				
				//allocate the array
				byteArray = buf.getBuffer(byteSize);
				//byteArray = new byte[byteSize];				
				
				//read that many bytes
				this.buffInput.read(byteArray, 0, byteSize);
			}
			catch(EOFException e)
			{
	  		this.buffInput.close();
	  		this.readClosed = true;			
	  		
	  		return null;
			}
		}
			
		//get as byte array stream
  	SequenceMinHashes seqHashes = SequenceMinHashes.fromByteStream(new DataInputStream(new ByteArrayInputStream(byteArray)), this.offset);

  	return seqHashes;
  }
}
