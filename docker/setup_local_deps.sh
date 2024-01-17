#!/bin/bash

mkdir /dependencies
cd /dependencies

# NIS deps
git clone https://github.com/CLIUtils/CLI11.git
git clone https://github.com/catchorg/Catch2.git
git clone https://gitlab.com/libeigen/eigen.git
git clone https://github.com/njoy/ENDFtk.git

# ENDFtk deps
git clone https://github.com/gabime/spdlog.git
cd spdlog
git checkout a51b4856377a71f81b6d74b9af459305c4c644f8
cd ..

git clone https://github.com/njoy/disco.git
cd disco
git checkout 2606933a854bb0269c4ec37143e1236797e27838
cd ..

git clone https://github.com/njoy/header-utilities.git
cd header-utilities
git checkout cc2610fee15e255c151e8e22aca1e8b3d1a96b39
cd ..

git clone https://github.com/njoy/Log.git
cd Log
git checkout 52962b7796afe37ef1d8f7edb4bf9ecb1b868d15
cd ..

git clone https://github.com/ericniebler/range-v3.git
cd range-v3
# v0.11.0
git checkout 4989f3e9ff2efee1852942bb9328ef121369ba02
cd ..