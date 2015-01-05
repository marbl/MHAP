/* 
 * MHAP package
 * 
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2014 by Konstantin Berlin and Sergey Koren
 * University Of Maryland
 * 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */
package edu.umd.marbl.mhap.general;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import edu.umd.marbl.mhap.utils.FastAlignRuntimeException;

public abstract class AbstractSequenceSearchMain<S extends AbstractMatchSearch<H>, H extends SequenceHashes>
{
	private final String processFile;
	private final String inFile;
	private final String toFile;
	private final boolean noSelf;
	protected final int numThreads;
	
	protected abstract AbstractSequenceHashStreamer<H> getSequenceHashStreamer(String file, int offset) throws IOException;
	protected abstract S getMatchSearch(AbstractSequenceHashStreamer<H> matchSearch) throws IOException;
	protected abstract void outputFinalStat(S matchSearch);
	
	public AbstractSequenceSearchMain(String processFile, String inFile, String toFile, boolean noSelf, int numThreads)
	{
		this.processFile = processFile;
		this.inFile = inFile;
		this.toFile = toFile;
		this.noSelf = noSelf;
		this.numThreads = numThreads;
	}
	
	public void computeMain() throws IOException
	{
		long startTotalTime = System.nanoTime();		
		long startTime = System.nanoTime();
		long processTime = System.nanoTime();
		
		//if processing a directory
		if (this.processFile!=null && !this.processFile.isEmpty())
		{
			File file = new File(this.processFile);			
			if (!file.exists())
				throw new FastAlignRuntimeException("Process file does not exist.");

			if (this.toFile==null || this.toFile.isEmpty())
				throw new FastAlignRuntimeException("Target directory must be defined.");
			
			File toDirectory = new File(this.toFile);
			if (!toDirectory.exists() || !toDirectory.isDirectory())
				throw new FastAlignRuntimeException("Target directory doesn't exit.");
			
			//allocate directory files
			ArrayList<File> processFiles = new ArrayList<>();
			
			//if not dictory just add the file
			if (!file.isDirectory())
			{
				processFiles.add(file);
			}
			else
			{			
				//read the directory content
				File[] fileList = file.listFiles(new FilenameFilter()
				{				
					@Override
					public boolean accept(File dir, String name)
					{
						if (!name.startsWith("."))
							return true;
						
						return false;
					}
				});
				
				for (File cf : fileList)
					processFiles.add(cf);				
			}
			
			for (File pf : processFiles)
			{
				startTime = System.nanoTime();
				
				AbstractSequenceHashStreamer<H> seqStreamer = getSequenceHashStreamer(pf.getAbsolutePath(), 0);
				
				String outputString = pf.getName();
				int i = outputString.lastIndexOf('.');
				if (i>0)
					outputString = outputString.substring(0, i);
				
				//combine with the directory name
				outputString = toDirectory.getPath()+File.separator+outputString+".dat";
				
				//store the file to disk
				seqStreamer.writeToBinary(outputString, false, this.numThreads);
				
				System.err.println("Processed "+seqStreamer.getNumberProcessed()+" sequences (fwd and rev).");
				System.err.println("Read, hashed, and stored file "+pf.getPath()+" to "+outputString+".");
				System.err.println("Time (s): " + (System.nanoTime() - startTime)*1.0e-9);
			}
			
			System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);

			return;
		}
		
		// read and index the kmers
		int seqNumberProcessed = 0;
				
		//create search object
		AbstractSequenceHashStreamer<H> seqStreamer = getSequenceHashStreamer(this.inFile, seqNumberProcessed);
		S hashSearch = getMatchSearch(seqStreamer);

		seqNumberProcessed += seqStreamer.getNumberProcessed()/2;
		System.err.println("Processed "+seqStreamer.getNumberProcessed()+" unique sequences (fwd and rev).");
		System.err.println("Time (s) to read and hash from file: " + (System.nanoTime() - processTime)*1.0e-9);

		long startTotalScoringTime = System.nanoTime();

		//System.err.println("Press Enter...");
		//System.in.read();
		
		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		if (this.toFile==null || this.toFile.isEmpty())
		{
			startTime = System.nanoTime();
			hashSearch.findMatches();
			System.err.println("Time (s) to score and output to self: " + (System.nanoTime() - startTime)*1.0e-9);
		}
		else
		{
			File file = new File(this.toFile);
			
			if (!file.exists())
				throw new FastAlignRuntimeException("To-file does not exist.");
			
			ArrayList<File> toFiles = new ArrayList<>();
			
			//if not dictory just add the file
			if (!file.isDirectory())
			{
				toFiles.add(file);
			}
			else
			{			
				//read the directory content
				File[] fileList = file.listFiles(new FilenameFilter()
				{				
					@Override
					public boolean accept(File dir, String name)
					{
						if (!name.startsWith("."))
							return true;
						
						return false;
					}
				});
				
				for (File cf : fileList)
					toFiles.add(cf);
			}

			//sort the files in alphabetical order
			Collections.sort(toFiles);

			//first perform to self
			startTime = System.nanoTime();
			if (!this.noSelf)
			{
				hashSearch.findMatches();
				System.out.flush();
				System.err.println("Time (s) to score and output to self: " + (System.nanoTime() - startTime)*1.0e-9);
			}

			//no do to all files
			for (File cf : toFiles)
			{			
				// read and index the kmers
				seqStreamer = getSequenceHashStreamer(cf.getAbsolutePath(), seqNumberProcessed);
				System.err.println("Opened fasta file "+cf.getCanonicalPath()+".");
	
				//match the file
				startTime = System.nanoTime();
				hashSearch.findMatches(seqStreamer);
				
				//flush to get the output
				System.out.flush();
				
				seqNumberProcessed += seqStreamer.getNumberProcessed();
				System.err.println("Processed "+seqStreamer.getNumberProcessed()+" to sequences.");
				System.err.println("Time (s) to score, hash to-file, and output: " + (System.nanoTime() - startTime)*1.0e-9);
			}
		}
		
		//flush output
		System.out.flush();
		
		//output time
		System.err.println("Total scoring time (s): " + (System.nanoTime() - startTotalScoringTime)*1.0e-9);
		System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);
		
		//output final stats
		outputFinalStat(hashSearch);
	}
}
