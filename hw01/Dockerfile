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

RUN ./scripts/run_linters.sh

RUN ./scripts/build.sh

RUN ./scripts/run_tests.sh

RUN ./scripts/generate_tests_cov_report.sh