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

## clone matplotlib-cpp
repository_matplotlib="https://github.com/lava/matplotlib-cpp.git"
localFolder_matplotlib="External/matplotlib"
git clone "$repository_matplotlib" "$localFolder_matplotlib"

## clone json
repository_json="https://github.com/nlohmann/json.git"
localFolder_json="External/Json"
git clone --depth 1 "$repository_json" "$localFolder_json"

## clone eigne
repository_eigen="https://gitlab.com/libeigen/eigen.git"
localFolder_eigen="External/Eigen"
git clone --depth 1 "$repository_eigen" "$localFolder_eigen"

