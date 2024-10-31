FROM registry.fedoraproject.org/fedora:36

ARG username=vscode

RUN yum update -y && \
    yum upgrade -y && \
    yum install \
        gcc \
        gfortran \
        g++ \
        gsl \
        gsl-devel \
        texlive \
        ImageMagick \
        dvipng \
        coreutils \
        qtchooser \
        cmake \
        git \
        ninja-build \
        eigen3-devel \
        tbb \
        tbb-devel \
        -y 

# Installs maquao (for profiling/development)
RUN cd /usr/bin && wget https://maqao.org/maqao_archive/maqao.intel64.2.17.0.tar.xz &&	\
    tar -xvf maqao.intel64.2.17.0.tar.xz && cp maqao.intel64.2.17.0/maqao.intel64 maqao && \
    rm -rf maqao.intel64.2.17.0.tar.xz maqao.intel64.2.17.0/

# Installs microbenchmark (for profiling/development)
RUN git clone https://gitlab.in2p3.fr/CTA-LAPP/PHOENIX_LIBS/MicroBenchmark.git && \
    cd MicroBenchmark && \
    ./install.sh /usr/ && \
    cd .. && \
    rm -fr MicroBenchmark

RUN useradd ${username}

USER ${username}

WORKDIR /home/${username}

RUN git clone https://github.com/docbrown1955/marty-public.git && \
    cd marty-public && \
    git checkout master && \
    mkdir bin include lib build 

RUN cd /home/${username}/marty-public/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=.install -G Ninja && \
    cmake --build . && \
    cmake --install .

ENV INSTALLMARTYPATH="/home/${username}/marty-public/build/.install"
ENV PATH=$INSTALLMARTYPATH/bin:$PATH
ENV LD_LIBRARY_PATH=$INSTALLMARTYPATH/lib
ENV LIBRARY_PATH=$INSTALLMARTYPATH/lib
ENV C_INCLUDE_PATH=$INSTALLMARTYPATH/include
ENV CPATH=$INSTALLMARTYPATH/include

# Installing DarkPACK with default options
COPY source.tar.gz /home/${username}/
RUN mkdir darkpack-cmake && \
    tar -xf source.tar.gz -C darkpack-cmake && \
    cd darkpack-cmake && \
    cmake -S . -B ./build \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=./build/.install \
        -G Ninja && \
    cmake --build ./build --target install_headers_all && \
    ninja -C build

CMD ["sleep", "infinity"]