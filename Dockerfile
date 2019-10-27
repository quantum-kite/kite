FROM ubuntu:18.10

ARG UNAME=user
ARG UID=1000
ARG GID=1000

RUN apt-get update
RUN apt-get -y install cmake g++
RUN apt-get -y install valgrind gdb strace
RUN apt-get -y install libhdf5-dev libeigen3-dev

RUN groupadd -g $GID -o $UNAME && useradd -m -u $UID -g $GID -o -s /bin/bash $UNAME

COPY . /home/${UNAME}/kite/
RUN chown -R ${UID}:${GID} /home/${UNAME}/kite/

USER $UNAME
WORKDIR /home/${UNAME}/kite/tools
RUN cmake . -DCMAKE_BUILD_TYPE=Debug && make
WORKDIR /home/${UNAME}/kite
RUN cmake . -DCMAKE_BUILD_TYPE=Debug && make
