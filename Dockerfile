FROM --platform=linux/amd64 rocker/geospatial:4.5

RUN apt-get update && apt-get upgrade -y
RUN R -e "install.packages(c('lidR', 'data.table', 'sf', 'jsonlite', 'argparse', 'rlas'))"

RUN mkdir -p /in && mkdir -p /out && mkdir -p /src/geodata
COPY src /src
COPY geodata /src/geodata

WORKDIR /src
CMD ["Rscript", "run.R"]