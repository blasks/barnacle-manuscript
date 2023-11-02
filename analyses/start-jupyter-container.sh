#!/bin/bash

# The purpose of this script is to start a jupyter notebook session using
# a software container.

pushd ../

docker run -v "${PWD}":/home -p 8888:8888 barnacle jupyter notebook --ip='*' --port=8888 --allow-root
