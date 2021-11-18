#!/bin/bash

COLLECTION="dao"
IMAGE="opencadc/${COLLECTION}2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run --rm --name ${COLLECTION}_run_vo --user $(id -u):$(id -g) -e HOME=/usr/src/app -v /daodomedata/archive/182:/data -v ${PWD}:/usr/src/app/ ${IMAGE} ${COLLECTION}_run_vo || exit $?

date
exit 0
