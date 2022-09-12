# SPPARKS grain growth simulations
This is a fork of the SPPARKS (11 Nov 2009) software package.
The official SPPARKS repository is located here: [https://github.com/spparks/spparks](https://github.com/spparks/spparks).

This project implements the SPPARKS portion of the *candidate grain* simulations as seen in DeCost & Holm [[1]](#1). The Python portion, which assigns grain types and provides some tools for post-processing simulation results, is in a separate repository: https://github.com/holmgroup/spparks-meso.





--------------------------------------------------------------------------

#  Installation 
## Docker
[Docker](https://www.docker.com/) is the preferred method of installation. The image can be built by simply running:
```bash
$ docker build -t spparks-candidate-grains .
```
## Singularity
Many HPC systems use [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) instead of Docker. The most straightforward way to build the singularity image is to build and export the image from a machine that has Docker, and then use singularity to convert it to the correct format. After building the image with the above command, export the container to a file:
```bash
$ docker save -o spparks-candidate-grains.tar spparks-candidate-grains:latest
```
Next, copy the image to a machine with singularity. The image can be converted to a singularity **.sif** file with the following command:

```bash
$ singularity build spparks-candidate-grains.sif docker-archive:./spparks-candidate-grains.tar
```

# Running a premade input deck.

SPPARKS provides several sample input decks in the `examples/` folder in this repository. The simple ising simulation runs very quickly and can be used to verify the installation works.
## Run with Docker

```bash
$ docker run --rm \
         -v $(pwd)/examples/ising/in.ising:/mnt/in.ising:ro \
         spparks-candidate-grains \
         spparks -in /mnt/in.ising > out.ising
```
  Summary of `docker run` arguments:
  - `--rm`: removes container after it executes, preventing clutter
  - `-v $(pwd)/examples/ising/in.ising:/mnt/in.ising:ro`  mounts the input file `in.ising` into the container as a read-only file. Note that the mount needs an absolute path.
  - `spparks-candidate-grains`: name of the docker image to run.
  - `spparks -in in.ising > results/out.ising`: Run spparks with the `in.ising` input deck. In this example, SPPARKS prints the results to the terminal. The results are saved to a file by redirecting the terminal output to `out.ising`.

Note: more advanced input decks will produce additional output files that cannot be handled through terminal redirection. In this case, you can either mount the output folder to the container (`-v path/on/host:path/in/container:rw`) or use [docker cp](https://docs.docker.com/engine/reference/commandline/cp/) to copy the outputs back to the host.

## Run with Singularity
Running with Singularity is very straightforward:
```bash
$ singularity run spparks-candidate-grains.sif spparks -in \
               examples/ising/in.ising > out.ising
```
Note that Singularity bind mounts `/home/$USER`, `/tmp`, and `$PWD` into the container at runtime, so no explicit mount commands are needed, unlike with Docker. After the container finishes executing, `out.ising` will be in the current directory.

# Building a new input deck
This repository is a fork of the 11 Nov 2009 release of [SPPARKS](https://github.com/spparks/spparks). The SPPARKS documentation covers the majority of the commands. This version provides some additional :

 - `app_style potts/ori`: For running candidate grains simulations with nonuniform boundary mobility from DeCost [[1]](#1).
 - `read_dream3d`: Read in in initial structure generated with [meso](https://github.com/holmgroup/spparks-meso).
 - `stats`: Save grain sizes at desired time interval in h5 file.
 - dream3d format for `dump`: Write results to dream3d file that can be processed with [meso](https://github.com/holmgroup/spparks-meso).

Sample input decks for candidate grains simulations are provided in the meso repository: https://github.com/holmgroup/spparks-meso. The meso repository also provides a command line interface for easily running candidate grains simulations.
# References 
<a id="1">[1]</a>
DeCost, B.L., Holm, E.A. Phenomenology of Abnormal Grain Growth in Systems with Nonuniform Grain Boundary Mobility. *Metall Mater Trans A* 48, 2771–2780 (2017).  [doi:10.1007/s11661-016-3673-6](https://doi.org/10.1007/s11661-016-3673-6)
