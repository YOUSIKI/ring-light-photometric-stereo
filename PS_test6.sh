#!/usr/bin/env bash

docker run --rm \
    --shm-size=1024M \
    -v "$PWD":/workspace \
    -w /workspace \
    -it \
    ring-light-photometric-stereo:latest \
    octave -q --eval "PS_test6" && echo "Execution completed." || echo "Execution failed."