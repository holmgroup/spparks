FROM ubuntu:18.04

RUN apt-get update && apt-get -y --no-install-recommends install \
	build-essential \
	libhdf5-dev \
	libeigen3-dev \
	&& apt-get clean \
	&& apt-get autoremove \
	&& rm -rf /var/lib/apt/lists*
	
# needed so linker can find right libraries to build spparks
RUN ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial.so /usr/lib/libhdf5.so
RUN ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so /usr/lib/libhdf5_hl.so

WORKDIR /home/spparks

COPY src src

RUN cd src && make stubs
RUN cd src && make h5

RUN ln -s /home/spparks/src/spk_h5 /usr/bin/spparks

RUN groupadd --gid=1000 spparks \
    && useradd --uid 1000 --gid 1000  spparks
