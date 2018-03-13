#!/bin/bash

# Change to directory of this build file
cd "$(dirname "$0")"

# Build base Docker image
docker build -t assemblyrepair:latest -f Dockerfile .
