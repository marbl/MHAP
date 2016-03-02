# MHAP

MinHash alignment process (MHAP pronounced MAP): locality sensitive hashing to detect overlaps and utilities. This is the development branch, please use the [latest tagged](https://github.com/marbl/MHAP/releases/tag/v2.0).

## Build

You must have a recent  [JDK](http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html "JDK") and [Apache Maven](http://maven.apache.org/ "MAVEN") available. To checkout and build run:

    git clone https://github.com/marbl/MHAP.git
    cd MHAP
    maven install
    
For a quick user-quide, run:

    cd target
    java -jar mhap-2.0.jar

## Docs
For the full documentation information please see http://mhap.readthedocs.org/en/

## License

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
