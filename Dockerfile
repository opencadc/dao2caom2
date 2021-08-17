FROM opencadc/matplotlib:3.9-slim

RUN apt-get update --no-install-recommends && apt-get dist-upgrade -y && \
    apt-get install -y build-essential \
                       git \
                       imagemagick && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN pip install cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    caom2utils \
    importlib-metadata \
    python-dateutil \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

ARG CAOM2_BRANCH=master
ARG CAOM2_REPO=opencadc
ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN git clone https://github.com/opencadc/cadctools.git && \
    cd cadctools && \
    pip install ./cadcdata && \
    cd ..

RUN git clone https://github.com/${CAOM2_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${CAOM2_BRANCH} && \
    pip install ./caom2utils && \
    cd ..


RUN pip install git+https://github.com/${PIPE_REPO}/caom2pipe@${PIPE_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/dao2caom2@${PIPE_BRANCH}#egg=dao2caom2

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
