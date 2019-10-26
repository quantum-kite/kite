FROM ubuntu:18.04

RUN apt-get update
RUN apt-get -y install cmake g++
RUN apt-get -y install valgrind gdb strace
RUN apt-get -y install libhdf5-dev libeigen3-dev

COPY . /root/kite/
WORKDIR /root/kite
RUN cmake . && make install
WORKDIR /root/kite/tools
RUN cmake . && make install
