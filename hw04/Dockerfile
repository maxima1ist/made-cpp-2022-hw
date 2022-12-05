FROM ubuntu:latest

RUN apt update && \
    apt install -y --no-install-recommends \
    gdb \
    git \
    build-essential \
    pkg-config \
    cmake \
    clang-format \
    libgtest-dev \
    valgrind \
    libpcsclite-dev \
    lcov \
    cppcheck \
    libboost-all-dev\
    aptitude \
    pip \
    unzip && \
    aptitude search boost && \
    apt clean && \
    rm -rf /var/lib/apt/lists

RUN cd /usr/src/googletest && \
    cmake . && \
    cmake --build . --target install

COPY . /hw04

WORKDIR /hw04

RUN pip install gdown && \
    cd /hw04 && \
    gdown https://drive.google.com/uc?id=1sahQsFxm_NYsCd4haShQDSkv8Men8VqG && \
    unzip -o data.zip && \
    rm data.zip

RUN ./scripts/build.sh

RUN ./scripts/run_tests.sh

RUN ./scripts/build.sh -p

RUN ./scripts/run_tests.sh