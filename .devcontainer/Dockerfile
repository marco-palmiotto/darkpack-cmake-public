FROM mcr.microsoft.com/devcontainers/cpp:1-debian-12

ARG username=vscode
ARG nthreads=8

# [Optional] Uncomment this section to install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libgsl-dev libgslcblas0 texlive texlive-luatex dvipng coreutils imagemagick \
    qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools \
    make cmake g++ gcc gfortran git neofetch

USER ${username}

WORKDIR /home/${username}

RUN git clone https://github.com/docbrown1955/marty-public.git && \
    cd marty-public && \
    git checkout master && \
    mkdir bin include lib build 

RUN cd /home/${username}/marty-public/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=.install && \
    make -j${nthreads} && \
    make install

# ARG INSTALLMARTYPATH_dockerfile="/home/${username}/marty-public/build"

# If $INSTALLMARTYPATH is not defined, we define it in the ~/.bashrc
# RUN cd /home/${username} && \
#     echo "# MARTY envirnmental variables" > .martyconfig && \
#     echo "export INSTALLMARTYPATH=${INSTALLMARTYPATH_dockerfile} # path where MARTY is installed">> .martyconfig && \
#     echo 'export PATH=$INSTALLMARTYPATH/bin:$PATH # Find the command grafed'>> .martyconfig && \
#     echo 'export LD_LIBRARY_PATH=$INSTALLMARTYPATH/lib:$LD_LIBRARY_PATH # Find libraries at run-time (Linux)'>>  .martyconfig && \
#     echo 'export LIBRARY_PATH=$INSTALLMARTYPATH/lib:$LIBRARY_PATH # Find libraries at compile-time' >> .martyconfig && \
#     echo 'export C_INCLUDE_PATH=$INSTALLMARTYPATH/include:$C_INCLUDE_PATH # Find headers at compile time' >> .martyconfig && \
#     echo 'export CPATH=$INSTALLMARTYPATH/include:$CPATH # Find headers at compile time' >> .martyconfig && \
#     cat .martyconfig >> .bashrc


ENV INSTALLMARTYPATH="/home/${username}/marty-public/build/.install"
ENV PATH=$INSTALLMARTYPATH/bin:$PATH
ENV LD_LIBRARY_PATH=$INSTALLMARTYPATH/lib
ENV LIBRARY_PATH=$INSTALLMARTYPATH/lib
ENV C_INCLUDE_PATH=$INSTALLMARTYPATH/include
ENV CPATH=$INSTALLMARTYPATH/include

CMD ["sleep", "infinity"]