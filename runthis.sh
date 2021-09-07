#!/bin/sh

# Setting up paths
echo 'export PATH="${HOME}/bin:$PATH"' >> ${HOME}/.zshrc
echo 'export PYTHONPATH="${HOME}:pythonlib:${PYTHONPATH}"' >> ${HOME}/.zshrc


# creating directories
mkdir ${HOME}/bin
mkdir ${HOME}/pythonlib


# copy files
cp -r ./lib/* ${HOME}/pythonlib/.
cp -r ./bin/* ${HOME}/bin/.
