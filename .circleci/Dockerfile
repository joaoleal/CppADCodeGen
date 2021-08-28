FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive TZ=Europe/Lisbon \ 
    apt-get install -y cmake \
                       cmake-curses-gui \
                       libclang-8-dev \
                       llvm-8-dev \
                       clang-8 \
                       gcc \
                       g++ \
                       wget \
                       doxygen \
                       graphviz \
                       git && \
    DEBIAN_FRONTEND=noninteractive TZ=Europe/Lisbon \
    apt-get install -y texlive-extra-utils \
                       texlive-extra-utils \
                       texlive-fonts-extra \
                       texlive-latex-recommended && \
    apt-get autoclean && \
    apt-get autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
