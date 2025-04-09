# moments++
## fast and accurate prediction of two-locus statistics

Original publication: https://www.biorxiv.org/content/10.1101/2025.02.19.639084v1

### Building

INSTALL.txt contains instructions about dependencies.
Once these are properly installed, moments++ can be built with the following:

To build with native architecture settings (`-march=native`)
```
mkdir build
cd build
cmake cmake -DNativeBuild=ON ..
make
make install
```

Alternatively, it is possible to use container images. 
These are straightforward to obtain but the executable is not optimized and should be slower to run.

To get the docker image for the main branch:
```
docker pull ghcr.io/gvbarroso/momentspp/momentspp:latest
```

To get the docker image for the release version
```
docker pull ghcr.io/gvbarroso/momentspp/momentspp:v0.1
```

To build a local apptainer (or singularity) image:
```
apptainer build momentspp.sif docker://ghcr.io/gvbarroso/momentspp/momentspp:v0.1
```

### Usage

After a successful build

```
momentspp param=opt.bpp
```

The apptainer can be run with 
```
./momentspp.sif param=opt.bpp
```

And docker can be run with
```
docker run ghcr.io/gvbarroso/momentspp/momentspp:latest param=opt.bpp
```

The 'doc' directory (and examples therein) contains instructions on how to write the options file and the associated yaml file which specifies the evolutionary model of interest.
