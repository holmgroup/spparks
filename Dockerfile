FROM debian:10.12-slim

# need g++, make, hdf5, and eigen to build spparks
# ln -s commands are needed for the linker to find the right libraries during the build
RUN apt-get update && apt-get -y --no-install-recommends install \
	build-essential \
	libhdf5-dev \
	libeigen3-dev && \
	apt-get clean && \
	apt-get autoremove && \
	rm -rf /var/lib/apt/lists* && \
	ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial.so /usr/lib/libhdf5.so && \
    ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so /usr/lib/libhdf5_hl.so

WORKDIR /home/spparks

# copy source code into container
COPY src src

# Stubs builds a dummy mpi library.
# Instead of running individual simulations in parallel with mpi,
# separate containers running different simulations can be run in parallel.
# Then build spparks (make h5) and put the binary in /usr/bin
RUN cd src && make stubs && make h5 && mv spk_h5 /usr/bin/spparks
