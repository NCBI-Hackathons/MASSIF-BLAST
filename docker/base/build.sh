#!/bin/bash

# Change to directory of this build file
cd "$(dirname "$0")"

# Build base Docker image
docker build -f Docker file -t assemblyrepair:latest
