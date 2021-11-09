#!/bin/bash

COLLECTION="dao"
IMAGE="opencadc/${COLLECTION}2caom2"

echo "Get image ${IMAGE}"
docker pull ${IMAGE} || exit $?

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Run image ${IMAGE}"
docker run --rm --name ${COLLECTION}_run  --user $(id -u):$(id -g) -e HOME=/usr/src/app -v /daodomedata/archive/182:/data -v ${PWD}:/usr/src/app/ ${IMAGE} ${COLLECTION}_run_state || exit $?

date
exit 0
