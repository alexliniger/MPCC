#!/bin/sh
## Copyright 2019 Alexander Liniger

## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at

##     http://www.apache.org/licenses/LICENSE-2.0

## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
###########################################################################
###########################################################################
## Install dependencies
set -e

## clone blasfeo
repository_blasfeo="https://github.com/giaf/blasfeo.git"
localFolder_blasfeo="External/blasfeo"
git clone "$repository_blasfeo" "$localFolder_blasfeo"
## clone hpipm
repository_hpipm="https://github.com/giaf/hpipm.git"
localFolder_hpipm="External/hpipm"
git clone "$repository_hpipm" "$localFolder_hpipm"
## clone matplotlib-cpp
repository_matplotlib="https://github.com/lava/matplotlib-cpp.git"
localFolder_matplotlib="External/matplotlib"
git clone "$repository_matplotlib" "$localFolder_matplotlib"
## clone eigne
repository_eigen="https://gitlab.com/libeigen/eigen.git"
localFolder_eigen="External/Eigen"
git clone "$repository_eigen" "$localFolder_eigen"
## clone json
repository_json="https://github.com/nlohmann/json.git"
localFolder_json="External/Json"
git clone "$repository_json" "$localFolder_json"


cd External/blasfeo
mkdir -p build
mkdir -p lib
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$(realpath ../lib)
make
make install

cd ../../hpipm
mkdir -p build
mkdir -p lib
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$(realpath ../lib) -DBLASFEO_PATH=$(realpath ../../blasfeo/lib)
make
make install
