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
 * Copyright (c) 2015 by Konstantin Berlin and Sergey Koren
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
package edu.umd.marbl.mhap.utils;

import java.util.LinkedList;
import java.util.ListIterator;

import edu.umd.marbl.mhap.impl.MhapRuntimeException;

public final class RandomSequenceGenerator
{
	private MersenneTwisterFast randGenerator;

	public RandomSequenceGenerator()
	{
		this.randGenerator = new MersenneTwisterFast();
	}

	public RandomSequenceGenerator(int seed)
	{
		this.randGenerator = new MersenneTwisterFast(seed);
	}

	private final char getRandomBase(Character toExclude)
	{
		Character result = null;

		while (result == null)
		{
			double base = this.randGenerator.nextDouble();
			if (base < 0.25)
			{
				result = 'A';
			}
			else if (base < 0.5)
			{
				result = 'C';
			}
			else if (base < 0.75)
			{
				result = 'G';
			}
			else
			{
				result = 'T';
			}

			if (toExclude != null && toExclude.equals(result))
			{
				result = null;
			}
		}

		return result;
	}
	
	public String generateRandomSequence(int length)
	{
		StringBuilder str = new StringBuilder(length);
		
		for (int iter=0; iter<length; iter++)
			str.append(getRandomBase(null));
		
		return str.toString();
	}
	
	//0.1188 0.0183 0.0129
	public String addPacBioError(String str)
	{
		return addError(str, 0.1188, 0.0183, 0.0129);
	}

	public String addError(String str, double insertionRate, double deletionRate, double substitutionRate)
	{
		if (insertionRate < 0.0 || deletionRate < 0.0 || substitutionRate < 0.0)
			throw new MhapRuntimeException("Error rate cannot be negative.");
		
		if (insertionRate+deletionRate+substitutionRate>1.00001)
			throw new MhapRuntimeException("Error rate must be less than or equal to 1.0.");

		double errorRate = insertionRate + deletionRate + substitutionRate;

		// use a linked list for insertions
		LinkedList<Character> modifiedSequence = new LinkedList<>();
		for (char a : str.toCharArray())
			modifiedSequence.add(a);

		// now mutate
		ListIterator<Character> iter = modifiedSequence.listIterator();
		while (iter.hasNext())
		{
			char i = iter.next();

			if (randGenerator.nextDouble() < errorRate)
			{
				double errorType = randGenerator.nextDouble();
				if (errorType < substitutionRate)
				{ // mismatch
					// switch base

					iter.set(getRandomBase(i));

					i++;
				}
				else if (errorType < insertionRate + substitutionRate)
				{ // insert

					iter.previous();
					iter.add(getRandomBase(null));
				}
				else
				{ // delete

					iter.remove();
				}
			}
			else
			{
				// i++;
			}
		}

		StringBuilder returnedString = new StringBuilder(modifiedSequence.size());
		for (char c : modifiedSequence)
			returnedString.append(c);

		return returnedString.toString();
	}

}
