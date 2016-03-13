# AdaptTools
This directory contains a set of tools for mesh adaptation intwo and three dimensions. They can be used and combined in a Linux/Unix/Posix-compliant shell script wrapper or called as functions from C/C++/F programs.

#### Mesh adaptation

#### Installation
1. Install the [ICS Commons Library](https://github.com/ICStoolbox/Commons) on your system. 
Please refer to the instructions provided on the ICS Commons Library page in order to install this library.

2. download the zip archive of AdaptTools or clone this repository:

   ` git clone https://github.com/ICStoolbox/AdaptTools.git `

   navigate to the downloaded directory: 

   ` cd AdaptTools `

   create a build directory and compile with cmake
   ```
   mkdir build
   cd build
   cmake ..
   make
   make install
   ```

#### Authors & contributors
* these tools have been developed and contributed over the years by several contributors. Main developers are Charles Dapogny (Université J. Fourier), Cécile Dobrzynski (Université de Bordeaux, INRIA) and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
AdaptTools is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

