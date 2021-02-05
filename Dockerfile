FROM opencadc/matplotlib:3.8-slim

RUN apt-get update -y && apt-get dist-upgrade -y && \
    apt-get install -y build-essential \
                       git \
                       imagemagick && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN pip install cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    ftputil \
    importlib-metadata \
    pytz \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN pip install git+https://github.com/${PIPE_REPO}/caom2pipe@${PIPE_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/dao2caom2@${PIPE_BRANCH}#egg=dao2caom2

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
