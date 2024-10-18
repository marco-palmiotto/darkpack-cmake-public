FROM mcr.microsoft.com/devcontainers/cpp:1-debian-12

ARG username=vscode

# Installing packages for MARTY's dependencies
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libgsl-dev libgslcblas0 texlive texlive-luatex dvipng coreutils imagemagick \
    qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools \
    make cmake g++ gcc gfortran git neofetch \
    autoconf automake libtool libgoogle-glog-dev libgflags-dev libeigen3-dev

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