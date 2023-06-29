# Use the official Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Set the geographic area environment variable
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Paris

# Set the working directory inside the container
WORKDIR /app

RUN apt update && \
    apt install -y build-essential wget && \
    apt install -y cmake && \
    apt install -y git && \
    apt install -y libgtest-dev

RUN apt-get update && apt-get install -y \
    python3 \
    python3-dev \
    python3-pip\
    python3-numpy

# Install Python packages
RUN pip3 install numpy==1.20.0 matplotlib

# Install Eigen library
RUN apt-get install -y libeigen3-dev

RUN git clone https://github.com/coin-or/CppAD.git
WORKDIR /app/CppAD
RUN mkdir build
WORKDIR /app/CppAD/build
RUN cmake ..
RUN make install

WORKDIR /app

RUN git clone https://github.com/joaoleal/CppADCodeGen.git
WORKDIR /app/CppADCodeGen
RUN mkdir build
WORKDIR /app/CppADCodeGen/build
RUN cmake ..
RUN make install

WORKDIR /app

RUN git clone https://github.com/giaf/blasfeo.git
WORKDIR /app/blasfeo
RUN mkdir build
WORKDIR /app/blasfeo/build
RUN cmake ..
RUN make install

WORKDIR /app
RUN git clone https://github.com/giaf/hpipm.git
WORKDIR /app/hpipm
RUN make static_library TARGET=GENERIC -j 4

# Create a directory inside the container to mount the external folder
RUN mkdir /app/MPCC

# Mount the external folder to /app/MPCC inside the container
VOLUME /app/MPCC

# Change to the MPCC directory
WORKDIR /app/MPCC
