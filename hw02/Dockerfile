FROM ubuntu:latest

RUN apt update && \
    apt install -y --no-install-recommends \
    gdb \
    git \
    build-essential \
    pkg-config \
    cmake \
    libgtest-dev \
    valgrind \
    libpcsclite-dev \
    lcov \
    cppcheck && \
    apt clean && \
    rm -rf /var/lib/apt/lists

RUN cd /usr/src/googletest && \
    cmake . && \
    cmake --build . --target install

COPY . /hw02

WORKDIR /hw02

RUN ./scripts/build.sh

RUN ./scripts/run_tests.sh