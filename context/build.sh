#!/bin/bash

docker build --no-cache -t "ghcr.io/imsb-uke/maxomod_pipeline-r-python:latest" .
docker build --no-cache -t "ghcr.io/imsb-uke/maxomod_pipeline-mofa:latest" --file  Dockerfile_mofa .
