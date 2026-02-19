**WHAT IS A DOCKERFILE?**
Dockerfiles are used to create Docker images.  Docker images are the 'recipes' used to build virtual environments called containers. These dockerfiles all use a base Linux image that has a particular compiler installed.  These base images were created by David Beckingsale and are created from this GitHub respository:
https://github.com/LLNL/radiuss-docker

**EXPLANATION OF DOCKERFILE CONTENTS**
The base images are accessed in these dockerfiles by the use of the 'FROM' command:
>FROM axom/compilers:&lt;compiler&gt;

The remaining dockerfile commands attempt to build the Axom third party libraries using the specified compiler. The following commands update the compiler and add Fortran compiler support:

>RUN sudo apt-get update -y

>RUN sudo apt-get install curl -fy
RUN sudo apt-get -qq install -y --no-install-recommends gfortran-7 && sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-7 100

>#(For Clang Builds - Need to restore stdc++ lib)

>RUN sudo apt-get install libstdc++-7-dev

Then the Axom and RAJA source is retrieved and RAJA is built and installed: (WORKDIR functions as a 'cd' command)
>RUN git clone --recursive https://github.com/LLNL/axom.git

>RUN git clone --recursive https://github.com/LLNL/RAJA.git

>WORKDIR "/home/axom/RAJA"

>RUN mkdir build

>WORKDIR "/home/axom/RAJA/build"

>RUN cmake ../. -DENABLE_TESTS=OFF

>RUN make 

>RUN sudo make install

Finally, the uberenv script is invoked to build the libraries and store them in the axom_tpls directory.  (The host config files produced during this process had been added to the Axom repo under "axom/host-configs/docker")

>sudo apt-get update -y

>WORKDIR "/home/axom/axom"

>RUN ./scripts/uberenv/uberenv.py --spack-config-dir=./scripts/uberenv/spack_configs/docker --spec=%&lt;compiler&gt;^mpich@3.0.4 --prefix=/home/axom/axom_tpls -k

Note: The mpich specification (^mpich@3.0.4) is necessary for clang builds.  Not specifying it currently results in the system attempting to use and patch version 3.2.1 of mpich, which fails. See https://github.com/spack/spack/pull/8320 and https://github.com/spack/spack/issues/8432
Also, mpich itself is used because openmpi's mpiexec will not run for the root user and causes multiple test failures.

Once completed, this container can be used to build and test the Axom application on a Docker-using CI web tool such as Azure Pipelines.

Already built Docker images containing these third party libraries can be accessed on Docker Hub at this location:
https://hub.docker.com/r/axom/tpls/tags

**BUILDING THE IMAGE**
Using a dockerfile to create an image requires a Docker client called Docker Desktop that provides a command line interface for Docker commands.  Once this client has been installed the image can be built in two ways.

1. BUILD INTERACTIVELY

   The first way is to run an interactive session of the base image followed by running the individual commands at the session prompt.

   Running the interactive session is accomplished with the "docker run -it <base image>" command. For example:

   >docker run -it axom/compilers:&lt;compiler&gt;

   This should start a container that provides a linux prompt.  The remaining commands can be entered in order (without the Docker 'RUN' command):

   >sudo apt-get update -y

   >sudo apt-get install curl -fy

   etc...

   Once all the commands have been completed the container can be exited and then converted into an image using the "docker commit" command.

   >docker commit &lt;container id&gt; &lt;image repository:tag&gt;

   The container id can be found with "docker container ls -a".
   The &lt;image repository:tag&gt; is chosen by the user to identify this image.

2. BUILD WITH DOCKERFILE

   The second way is to have the Docker engine build the image directly from the dockerfile.  This is done by placing the dockerfile in an empty directory and using the "docker build" command.

   >docker build -t &lt;image repository:tag&gt; &lt;dir with dockerfile&gt;

   The -t &lt;image repositag&gt; opt sets the repository:tag identifier of the image.

**PUSHING THE IMAGE TO DOCKER HUB**
Once created the image can be uploaded to Docker Hub with the "docker push" command.
>docker push &lt;image repository:tag&gt;
