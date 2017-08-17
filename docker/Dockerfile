FROM debian:jessie
MAINTAINER Markus Lux

RUN apt-get update
RUN apt-get -y install build-essential cmake git libboost-filesystem-dev libboost-program-options-dev libboost-system-dev hmmer2 libxml-simple-perl

# install acdc

RUN git clone https://github.com/mlux86/acdc.git /tmp/acdc
RUN mkdir /tmp/acdc/build
WORKDIR /tmp/acdc/build

RUN cmake ..
RUN make -j $(nproc)
RUN make install

RUN rm -rf /tmp/acdc

# install kraken

RUN git clone https://github.com/DerrickWood/kraken.git /tmp/kraken
WORKDIR /tmp/kraken
RUN ./install_kraken.sh /usr/local/bin
RUN rm -rf /tmp/kraken

# install rnammer

RUN mkdir /tmp/rnammer
ADD rnammer-1.2.src.tar.Z /tmp/rnammer
ADD rnammer.patch /tmp/rnammer
WORKDIR /tmp/rnammer

RUN tar xf rnammer-1.2.src.tar.Z
RUN rm rnammer-1.2.src.tar.Z
RUN patch < rnammer.patch
RUN mkdir /usr/local/share/rnammer
RUN cp -r * /usr/local/share/rnammer
RUN ln -s /usr/local/share/rnammer/rnammer /usr/local/bin/rnammer
RUN rm -rf /tmp/rnammer
WORKDIR /

# provide data shares and setup nginx compability

RUN mkdir /acdc

# start acdc per default

ENTRYPOINT ["/usr/local/bin/acdc", "-o", "/acdc", "-K", "/krakendb"]

# Run acdc
# docker run --name acdc -v /path/to/assemblies:/assemblies -v /path/to/minikraken_db:/krakendb acdc -i /assemblies/test.fasta

# Copy results from container to current working directory
# docker cp acdc:/acdc ./results