FROM opencadc/matplotlib:3.8-slim

RUN apt-get update
RUN apt-get install -y \
    build-essential \
    git

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

RUN apt-get install -y imagemagick

WORKDIR /usr/src/app

ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2pipe.git && \
    pip install ./caom2pipe
  
RUN git clone https://github.com/${OPENCADC_REPO}/dao2caom2.git && \
  cp ./dao2caom2/scripts/config.yml / && \
  cp ./dao2caom2/scripts/docker-entrypoint.sh / && \
  pip install ./dao2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]
