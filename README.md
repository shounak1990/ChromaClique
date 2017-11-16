# ChromaClique

ChromaClique finds nucleosome profiles from NOMe-sequencing reads
We present a novel method to reconstruct the structure of a viral quasispecies from NGS data.
Our approach can be used to:
 - reconstruct local error-corrected haplotypes and estimate their abundance
 - assemble full-length viral haplotypes
 - detect large deletions and insertions from paired-end data.


### INSTALL
Dependencies
ChromaClique depends on [boost](http://www.boost.org/), [gnu parallel](http://www.gnu.org/software/parallel/), and [cmake](http://www.cmake.org/). You can install them with a package manager of your choice.

Ubuntu:  
```
apt-get install libncurses5-dev cmake libboost-all-dev git build-essential zlib1g-dev parallel
```

OSX 10.8.x:
Please XCode and its command line tools, and with [macports](http://www.macports.org/):
```
port install cmake boost parallel
```

Windows:
ChromaClique has not been tested on Windows. The scripts depend on the bash shell, awk, and sed.  

Installation routine:
If you want to install ChromaClique to a non-standard directory, change it with `cmake -DCMAKE_INSTALL_PREFIX=<prefix-path> ..`
```bash
git clone https://github.com/shounak1990/ChromaClique
cd ChromaClique
git submodule init && git submodule update
mkdir build
cd build
cmake ..
make
make install
```
### Manual
Please visit https://shounak1990.github.io/chromaclique/manual.html for the manual.

### Contributions
 - [Armin Töpfer](http://www.armintoepfer.com)  
 - [Tobias Marschall](https://bioinf.mpi-inf.mpg.de/homepage/index.php?&account=marschal)
 - Bernhard Lang
 - Marcel Meyerheim

Contact
```
Shounak Chakraborty
shounak.chakraborty1990 (at) gmail.com
```
