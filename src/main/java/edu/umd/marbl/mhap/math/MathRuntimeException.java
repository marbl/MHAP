/* 
 * ARMOR package
 * 
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2012 by Konstantin Berlin 
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
package edu.umd.marbl.mhap.math;

/**
 * The Class MathException.
 */
public class MathRuntimeException extends RuntimeException
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 6939427297903213601L;

	/**
	 * Instantiates a new math exception.
	 */
	public MathRuntimeException()
	{
		super();
	}

	/**
	 * Instantiates a new math exception.
	 * 
	 * @param arg0
	 *            the arg0
	 */
	public MathRuntimeException(String arg0)
	{
		super(arg0);
	}

	/**
	 * Instantiates a new math exception.
	 * 
	 * @param arg0
	 *            the arg0
	 * @param arg1
	 *            the arg1
	 */
	public MathRuntimeException(String arg0, Throwable arg1)
	{
		super(arg0, arg1);
	}

	/**
	 * Instantiates a new math exception.
	 * 
	 * @param arg0
	 *            the arg0
	 */
	public MathRuntimeException(Throwable arg0)
	{
		super(arg0);
	}

}
