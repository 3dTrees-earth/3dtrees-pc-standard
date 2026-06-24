FROM --platform=linux/amd64 rocker/geospatial:4.5
SHELL ["/bin/bash", "-c"]
ARG RLAS_COMMIT=82cbba42f158d1dfc91efda3207923260a052564
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl patch && \
    rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages(c('lidR', 'data.table', 'sf', 'jsonlite', 'argparse', 'Rcpp'))"
COPY patches/rlas-read-all-extrabytes.patch /tmp/rlas-read-all-extrabytes.patch
RUN curl -L -sS "https://github.com/r-lidar/rlas/archive/${RLAS_COMMIT}.tar.gz" -o /tmp/rlas.tar.gz && \
    mkdir -p /tmp/rlas-src && \
    tar -xzf /tmp/rlas.tar.gz -C /tmp/rlas-src --strip-components=1 && \
    patch -d /tmp/rlas-src -p1 < /tmp/rlas-read-all-extrabytes.patch && \
    R CMD INSTALL /tmp/rlas-src && \
    rm -rf /tmp/rlas-src /tmp/rlas.tar.gz /tmp/rlas-read-all-extrabytes.patch

RUN mkdir -p /in && mkdir -p /out && mkdir -p /src/geodata
COPY src /src
COPY geodata /src/geodata

RUN chmod -R 777 /in /out /src

WORKDIR /src
CMD ["Rscript", "run.R"]
