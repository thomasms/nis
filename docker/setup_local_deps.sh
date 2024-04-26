#!/bin/bash

mkdir ./dependencies
cd ./dependencies

# NIS deps
git clone https://github.com/CLIUtils/CLI11.git
git clone https://github.com/catchorg/Catch2.git
git clone https://gitlab.com/libeigen/eigen.git

# ENDFtk - v1.0.0
git clone https://github.com/njoy/ENDFtk.git
cd ENDFtk
git checkout v1.0.0
cd ..

# ENDFtk deps
# tag: v1.11.0
git clone https://github.com/gabime/spdlog.git
cd spdlog
git checkout ad0e89cbfb4d0c1ce4d097e134eb7be67baebb36
cd ..

git clone https://github.com/njoy/disco.git
cd disco
git checkout 2606933a854bb0269c4ec37143e1236797e27838
cd ..

# tag: v0.2.0
git clone https://github.com/njoy/tools.git
cd tools
git checkout 25c9273d05601a9644feea6d7539250bf1d1c0dc
cd ..

# tag: 0.11.0
git clone https://github.com/ericniebler/range-v3.git
cd range-v3
# v0.11.0
git checkout 4989f3e9ff2efee1852942bb9328ef121369ba02
cd ..