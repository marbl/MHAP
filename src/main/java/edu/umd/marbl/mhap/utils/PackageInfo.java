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
package edu.umd.marbl.mhap.utils;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * The Class PackageInfo.
 */
public final class PackageInfo
{

	private static Properties properties = getProjectProperties();

	private static Properties getProjectProperties()
	{
		Properties initialProperties = new Properties();

		initialProperties = getFileProperties("/properties/mhap.properties", new Properties());
		Properties properties = getFileProperties("/mhap.properties", initialProperties);

		return properties;
	}

	private static Properties getFileProperties(String file, Properties originalProperties)
	{
		Properties property = new Properties(originalProperties);
		try
		{
			InputStream in = PackageInfo.class.getClass().getResourceAsStream(file);

			if (in != null)
			{
				property.load(in);
				in.close();
			}
		}
		catch (IOException e)
		{
			return property;
		}

		return property;
	}

	public static final String VERSION = properties.getProperty("component.version", "unknown");
	public static final String BUILD_TIME = properties.getProperty("buildtime", "unknown");

	public static Properties getProperties()
	{
		return properties;
	}

	public static String getVersionTag()
	{
		return VERSION + " " + BUILD_TIME;
	}
}
