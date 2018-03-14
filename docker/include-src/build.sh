#!/bin/bash

# Change to directory of this build file
cd "$(dirname "$0")"

cp -R ../../src/ src/

# Build base Docker image
docker build -t assemblyscripts:latest -f Dockerfile .
