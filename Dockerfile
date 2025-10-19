FROM --platform=linux/amd64 rocker/geospatial:4.5
SHELL ["/bin/bash", "-c"]
RUN apt-get update && apt-get upgrade -y
RUN R -e "install.packages(c('lidR', 'data.table', 'sf', 'jsonlite', 'argparse', 'rlas', 'rcpp'))"

RUN mkdir -p /in && mkdir -p /out && mkdir -p /src/geodata
COPY src /src
COPY geodata /src/geodata

RUN chmod -R 777 /in /out /src

WORKDIR /src
CMD ["Rscript", "run.R"]