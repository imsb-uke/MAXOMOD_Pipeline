#!/bin/bash

# 64GB = 68719476736 bytes

MX_DOCKERIMAGE="ghcr.io/imsb-uke/maxomod_pipeline-r-python"
MX_TAG="latest"

which docker

if [ "$?" = "0" ]; then
  docker pull "${MX_DOCKERIMAGE}:${MX_TAG}"
  docker run \
    --rm \
    --name "$USER-maxomod-pipeline-$$" \
    --volume "$PWD:/workdir:rw" \
    --workdir "/workdir" \
    --user `id -u`:`id -g` \
    --cpus 8 \
    --memory 68719476736 \
    "${MX_DOCKERIMAGE}:${MX_TAG}" \
    $@
else
  echo "docker not available, running in shell"
  $@
fi
