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

COPY . /hw01

WORKDIR /hw01

RUN ./run_linters.sh

RUN ./build.sh

RUN ./run_tests.sh

RUN ./generate_tests_cov_report.sh