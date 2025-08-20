#!/bin/bash

mkdir ./dependencies
cd ./dependencies

# NIS deps
git clone https://github.com/CLIUtils/CLI11.git
git clone https://github.com/catchorg/Catch2.git
git clone https://gitlab.com/libeigen/eigen.git

# ENDFtk - v1.1.0
git clone https://github.com/njoy/ENDFtk.git
cd ENDFtk
git checkout v1.1.0
cd ..

# ENDFtk deps
# tag: v1.11.0
git clone https://github.com/gabime/spdlog.git
cd spdlog
git checkout ad0e89cbfb4d0c1ce4d097e134eb7be67baebb36
cd ..

git clone https://github.com/fastfloat/fast_float
cd fast_float
git checkout f476bc713fda06fbd34dc621b466745a574b3d4c
cd ..

# tag: v0.3.0
git clone https://github.com/njoy/tools.git
cd tools
git checkout 368dbd9bd44754de616c46ffed0f80d2d16d8360
cd ..

# tag: 0.11.0
git clone https://github.com/ericniebler/range-v3.git
cd range-v3
# v0.11.0
git checkout 4989f3e9ff2efee1852942bb9328ef121369ba02
cd ..