#!/usr/bin/env bash

# get sicor project
rm -rf context/sicor
git clone git@git.gfz-potsdam.de:EnMAP/sicor.git ./context/sicor
# git clone git@git.gfz-potsdam.de:EnMAP/sicor.git --branch feature/improve_enmap --single-branch ./context/sicor
cd ./context/sicor || exit
git lfs pull
pip wheel --no-deps .
cd ../..
