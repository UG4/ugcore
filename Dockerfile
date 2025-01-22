# syntax=docker/dockerfile:1
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

# Install ug4 dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    git \
    cmake \
    build-essential

# Setup ughub
WORKDIR /opt/
RUN git clone https://github.com/UG4/ughub
ENV PATH="$PATH:/opt/ughub"

# Setup minimal ug4
WORKDIR /opt/ug4
ENV UG4_ROOT=/opt/ug4/
RUN ughub init
RUN ughub install ugcore

# Build ug4
WORKDIR ${UG4_ROOT}/build
RUN cmake .. -DENABLE_ALL_PLUGINS=OFF -DCOMPILE_INFO=OFF -DCMAKE_BUILD_TYPE=Release 
RUN make -j
RUN make install
