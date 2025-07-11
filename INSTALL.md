## Dependencies:
- A C++ compiler, both gcc and clang should work. Recent versions are required, with support for C++20
- CMake >= 3.1 for building
- The Boost C++ libraries >= 1.71, available in https://www.boost.org/
- The Threads library (usually provided along with the C++ compiler installation)
- The Bio++ libraries >= 3.0 (bpp-core, bpp-seq, bpp-phyl), available in https://github.com/BioPP
- The Eigen3 libraries, available at https://eigen.tuxfamily.org/dox/GettingStarted.html
- The yaml-cpp library (https://github.com/jbeder/yaml-cpp)

## Installation using an external build:

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local .. (for a local install)
make
make install
```

The 'momentspp' executable should then be copied to $HOME/.local/bin 

NOTE: we can request 

```
-DVERBOSE=true ..
```

for compiling in 'verbose' mode, where diagnostic files will be output to the same directory. Furthermore, in Linux machines

```
-DNativeBuiild=ON
```

will optimize compilation by exploiting knowledge of the machine's architecture (this can decrease runtime by 20%).

We can put it all together, e.g.

```
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local -DNativeBuild=ON -DVERBOSE=true ..
```

On a standard Linux workstation, Moments++ takes < 2 minutes to compile and install once all dependencies are met. Note that static linkage is not supported at the moment, but out-of-the-box executable versions can be obtained with container images (see README in the home directory).

## NOTE ABOUT COMPUTING CLUSTERS:

It is possible that the default versions of gcc and cmake on computing clusters will not meet the aforementioned requirements.
To fix this, one can execute, e.g.,

```
module load gcc/11.2.0
module load cmake/3.7.2
```

in the interactive session that is used to build Moments++ on the cluster.

I have also noticed that cmake may have problems finding the Boost installation when working on a cluster.
To work around this problem one can execute (from the external build directory):

```
ccmake ..
```

then press "t" to enter "advanced" mode and edit the BOOST lines at the top of the GUI by pressing "enter" and enter the correct path.

then press "c" to configure and "g" to generate the make file. If the make file is generated successfully, one can then run "make install".

Once the Boost libraries are eventually found by cmake, compiling Moments++ on a cluster may throw several warnings regarding Boost, but so far I have been able to ignore them without problems.

