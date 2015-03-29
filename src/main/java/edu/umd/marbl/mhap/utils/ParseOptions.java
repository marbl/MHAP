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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * The Class ParseOptions.
 */
public class ParseOptions
{
	public class Option<T extends Object>
	{
		protected final Class<T> objectClass;

		private final String description;

		private final String flag;

		private Object value;

		private boolean wasSet;

		private boolean required;

		@SuppressWarnings("unchecked")
		private Option(String flag, String description, T defaultValue)
		{
			this.flag = flag;
			this.description = description;
			this.value = defaultValue;
			this.wasSet = false;
			this.required = false;
			this.objectClass = (Class<T>) defaultValue.getClass();
		}

		private Option(String flag, String description, boolean required, Class<T> objectClass)
		{
			this.flag = parseFlag(flag);
			this.description = description;
			this.value = null;
			this.wasSet = false;
			this.required = required;
			this.objectClass = objectClass;
		}

		private String parseFlag(String flag)
		{
			flag = flag.trim();
			if (flag.startsWith("-"))
				flag = flag.substring(1);

			return flag;
		}

		public boolean isBoolean()
		{
			return Boolean.class.isAssignableFrom(this.objectClass);
		}

		public boolean isDouble()
		{
			return Double.class.isAssignableFrom(this.objectClass);
		}

		public boolean isInteger()
		{
			return Integer.class.isAssignableFrom(this.objectClass);
		}

		public boolean isString()
		{
			return String.class.isAssignableFrom(this.objectClass);
		}

		public boolean getBoolean()
		{
			return (Boolean) this.value;
		}

		public String getString()
		{
			if (this.value == null)
				return "";

			return (String) this.value;
		}

		public double getDouble()
		{
			return (Double) this.value;
		}

		public int getInteger()
		{
			return (Integer) this.value;
		}

		public void setValue(Object value)
		{
			if (!this.objectClass.isInstance(value))
				throw new RuntimeException("Incompatable value type for flag \"" + this.flag + "\" of type "
						+ this.objectClass.getName() + ": " + value.getClass().getName());

			this.value = value;
			this.wasSet = true;
		}

		public String getFlag()
		{
			return this.flag;
		}

		public boolean isSet()
		{
			return this.wasSet;
		}

		public boolean isRequired()
		{
			return this.required;
		}

		public Class<T> getType()
		{
			return this.objectClass;
		}
	}

	private final ArrayList<String> startText;
	private final HashMap<String, Option<?>> optionsByFlag;
	private final HashMap<String, Option<?>> optionsByName;

	public ParseOptions()
	{
		this.startText = new ArrayList<>();
		this.optionsByFlag = new HashMap<String, Option<?>>();
		this.optionsByName = new HashMap<String, Option<?>>();
		addOption("-h", "Displays the help menu.", false);
		addOption("--help", "Displays the help menu.", false);
		addOption("--version", "Displays the version and build time.", false);
	}
	
	public void addStartTextLine(String text)
	{
		this.startText.add(text);
	}

	public <T> void addOption(String flag, String description, T defaultValue)
	{
		flag = parseFlag(flag);
		if (this.optionsByFlag.get(flag) != null)
			return;

		Option<T> option = new Option<T>(flag, description, defaultValue);
		this.optionsByFlag.put(flag, option);
		this.optionsByName.put(flag, option);
	}
	
	public <T> void setOptions(String flag, T defaultValue)
	{
		Option<T> option = new Option<T>(flag, getFlag(flag).description, defaultValue);
		this.optionsByFlag.put(flag, option);
		this.optionsByName.put(flag, option);		
	}

	public <T> void addRequiredOption(String flag, String description, Class<T> objectClass)
	{
		flag = parseFlag(flag);

		Option<T> option = new Option<T>(flag, description, true, objectClass);
		this.optionsByFlag.put(flag, option);
		this.optionsByName.put(flag, option);
	}

	public String parseFlag(String flag)
	{
		return flag;
	}

	public boolean process(String[] args)
	{
		try
		{
			parse(args);
			if (needsHelp())
			{
				System.out.println(helpMenuString());
				return false;
			}
			else
			if (needsVersion())
			{
				System.out.println("MHAP Version = "+PackageInfo.VERSION+", Build time = "+PackageInfo.BUILD_TIME);
				return false;
			}

			checkParameters();
		}
		catch (Exception e)
		{
			System.out.println(e.getMessage());
			System.out.println(helpMenuString());
			return false;
		}

		return true;
	}

	public Option<?> get(String name) throws RuntimeException
	{
		Option<?> option = this.optionsByName.get(name);

		if (option == null)
			throw new RuntimeException("Invalid option name \"" + name + "\".");

		return option;
	}

	public Option<?> getFlag(String flag) throws RuntimeException
	{
		Option<?> option = this.optionsByFlag.get(flag);

		if (option == null)
			throw new RuntimeException("Invalid flag \"" + flag + "\".");

		return option;
	}

	@Override
	public String toString()
	{
		StringBuilder menuString = new StringBuilder();

		// sort the list
		ArrayList<String> list = new ArrayList<String>(this.optionsByFlag.keySet());
		Collections.sort(list);

		for (String key : list)
		{
			Option<?> currOption = this.optionsByFlag.get(key);
			menuString.append("" + currOption.flag + " = ");
			menuString.append("" + currOption.value);
			menuString.append("\n");
		}

		return menuString.toString();
	}

	public String helpMenuString()
	{
		StringBuilder menuString = new StringBuilder();
		for (String str : this.startText)
			menuString.append(str+"\n");

		// sort the list
		ArrayList<String> list = new ArrayList<String>(this.optionsByFlag.keySet());
		Collections.sort(list);

		for (String key : list)
		{
			Option<?> currOption = this.optionsByFlag.get(key);
			menuString.append("\t\t" + currOption.flag + ", ");
			if (currOption.isRequired())
				menuString.append("*required, ");
			else
			{
				if (currOption.isString())
					menuString.append("default = \"" + currOption.value + "\"");
				else
					menuString.append("default = " + currOption.value);
			}
			menuString.append("\n");
			menuString.append("\t\t\t" + currOption.description);
			menuString.append("\n");
		}

		return menuString.toString();
	}

	public void checkParameters()
	{
		for (Option<?> option : this.optionsByFlag.values())
			if (option.required && !option.wasSet)
				throw new RuntimeException("Required option flag \"" + option.flag + "\" was not set.");
	}

	public boolean needsHelp()
	{
		return get("--help").getBoolean() || get("-h").getBoolean();
	}
	
	public boolean needsVersion()
	{
		return get("--version").getBoolean();
	}


	public void parse(String[] args) throws RuntimeException
	{
		for (int iter = 0; iter < args.length; iter++)
		{
			String flag = args[iter].trim();

			if (!flag.startsWith("-"))
				throw new RuntimeException("Unknown parameter in command line: " + flag);

			flag = parseFlag(flag);

			Option<?> option = getFlag(flag);
			if (option == null)
				throw new RuntimeException("Unknown flag \"" + flag + "\".");
			else if (option.isBoolean())
				option.setValue(true);
			else if (iter + 1 < args.length && !args[iter + 1].startsWith("-"))
			{
				if (option.isDouble())
				{
					option.setValue(new Double(args[iter + 1]));
					iter++;
				}
				else if (option.isInteger())
				{
					option.setValue(new Integer(args[iter + 1]));
					iter++;
				}
				else if (option.isString())
				{
					option.setValue(args[iter + 1]);
					iter++;
				}
				else
					throw new RuntimeException("Cannot parse flag \"" + option.getFlag() + "\" of type "
							+ option.getType().getName() + ".");
			}
			else
				throw new RuntimeException("Not value provided for flag \"" + option.getFlag() + "\" of type "
						+ option.getType().getName() + ".");
		}
	}
}