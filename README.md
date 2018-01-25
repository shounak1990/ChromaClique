# ChromaClique

ChromaClique is a tool for the reconstruction of nucleosome profiles of individual or subpopulations of cells. Data from NOMe-sequencing  is used by ChromaClique for deconvolution of different nucleosome profiles (chromatypes) from cell subpopulations of one NOME-seq measurement. ChromaClique uses a maximal clique enumeration algorithm on a newly defined NOMe read graph that is able to group reads according to their nucleosome profiles.

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
 - [Marcel Schulz](https://bioinf.mpi-inf.mpg.de/homepage/index.php?&account=mschulz)  
 - [Tobias Marschall](https://bioinf.mpi-inf.mpg.de/homepage/index.php?&account=marschal)
 - [Stefan Canzar](http://www.genzentrum.uni-muenchen.de/research-groups/canzar/group-members/canzar-stefan/index.html)
 - [Shounak Chakraborty](http://www.genzentrum.uni-muenchen.de/research-groups/canzar/group-members/chakraborty-shounak/index.html)

Contact
```
Shounak Chakraborty
shounak.chakraborty1990 (at) gmail.com
```
