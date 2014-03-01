package com.secret.fastalign.general;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Locale;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Utils;

public class FastaData
{
	private final BufferedReader fileReader;
	private String lastLine;
	private AtomicLong numberProcessed;
	private boolean readFullFile;
	// length of sequences loaded
	private final ConcurrentLinkedQueue<Sequence> sequenceList;

	private static final String[] fastaSuffix = { "fna", "contigs", "contig", "final", "fasta", "fa" };

	private FastaData(ConcurrentLinkedQueue<Sequence> seqList)
	{
		this.sequenceList = new ConcurrentLinkedQueue<Sequence>(seqList);
		this.fileReader = null;
		this.lastLine = null;
		this.readFullFile = true;
		this.numberProcessed = new AtomicLong(this.sequenceList.size());
	}

	public FastaData(String file) throws IOException
	{
		try
		{
			this.fileReader = Utils.getFile(file, fastaSuffix);
		}
		catch (Exception e)
		{
			throw new FastAlignRuntimeException(e);
		}

		this.lastLine = null;
		this.readFullFile = false;
		this.numberProcessed = new AtomicLong(0);
		this.sequenceList = new ConcurrentLinkedQueue<Sequence>();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public synchronized FastaData clone()
	{
		// enqueue all the data
		try
		{
			enqueueFullFile();
		}
		catch (IOException e)
		{
			throw new FastAlignRuntimeException(e);
		}

		return new FastaData(this.sequenceList);
	}

	public synchronized Sequence dequeue() throws IOException
	{
		if (this.sequenceList.isEmpty())
			enqueueNextSequenceInFile();

		return this.sequenceList.poll();
	}

	public synchronized void enqueue(Sequence seq)
	{
		this.sequenceList.add(seq);
		this.numberProcessed.getAndIncrement();
	}

	public synchronized void enqueueFullFile() throws IOException
	{
		while (enqueueNextSequenceInFile())	{}
	}

	private synchronized boolean enqueueNextSequenceInFile() throws IOException
	{
		if (this.readFullFile)
			return false;

		// try to read the next line
		if (this.lastLine == null)
		{
			this.lastLine = this.fileReader.readLine();

			// there is no next line
			if (this.lastLine == null)
			{
				this.fileReader.close();
				this.readFullFile = true;
				return false;
			}
		}

		// process the header
		if (!this.lastLine.startsWith(">"))
			throw new FastAlignRuntimeException("Next sequence does not start with >. Invalid format.");

		// process the current header
		// parse the new header
		// header = this.lastLine.substring(1).split("[\\s]+", 2)[0];
		String header = "";
		this.lastLine = this.fileReader.readLine();

		StringBuilder fastaSeq = new StringBuilder();
		while (true)
		{
			if (this.lastLine == null || this.lastLine.startsWith(">"))
			{
				// enqueue sequence
				enqueue(new Sequence(fastaSeq.toString().toUpperCase(Locale.ENGLISH), new SequenceId(header)));

				if (this.lastLine == null)
				{
					this.fileReader.close();
					this.readFullFile = true;
				}

				return true;
			}

			// append the last line
			fastaSeq.append(this.lastLine);
			this.lastLine = this.fileReader.readLine();
		}

	}

	public long getNumberProcessed()
	{
		return this.numberProcessed.get();
	}

	public Sequence getSequence(SequenceId id)
	{
		if (id.isForward())
		{
			for (Sequence seq : this.sequenceList)
				if (seq.getId().equals(id))
					return seq;
		}

		id = id.complimentId();
		for (Sequence seq : this.sequenceList)
			if (seq.getId().equals(id))
				return seq.getReverseCompliment();

		return null;
	}

        public Sequence[] toArray() {
           return this.sequenceList.toArray(new Sequence[(int)this.getNumberProcessed()]);
        }

	public boolean isEmpty()
	{
		return this.sequenceList.isEmpty() && this.readFullFile;
	}
}
