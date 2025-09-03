FROM ubuntu:noble

ARG username=vscode

# Installing packages for MARTY's dependencies
RUN apt-get update && \
    apt-get -y install --no-install-recommends \
        coreutils \
        build-essential \
        gcc-14 \
        g++-14 \
        gfortran-14 \
        libgsl-dev \
        libgslcblas0 \
        qtbase5-dev \
        qtchooser \
        qt5-qmake \
        qtbase5-dev-tools \
        dvipng \
        lcov \
        cmake \
        texlive \
        texlive-latex-extra \
        texlive-luatex \
        imagemagick \
        ninja-build \
        git \
        clang-format && \
    apt-get clean

RUN apt-get install -y ca-certificates && \
    update-ca-certificates

RUN useradd -s /bin/bash ${username}

USER ${username}

WORKDIR /home/${username}

RUN git clone https://github.com/docbrown1955/marty-public.git && \
    cd marty-public && \
    git checkout separating_global && \
    mkdir bin include lib build 

ENV CC="/usr/bin/gcc-14"
ENV CXX="/usr/bin/g++-14"
ENV FC="/usr/bin/gfortran-14"

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

CMD ["sleep", "infinity"]
