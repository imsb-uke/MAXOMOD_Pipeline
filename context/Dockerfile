FROM rocker/verse:4.1.1

ENV POETRY_VIRTUALENVS_CREATE=true \
    POETRY_INSTALLER_PARALLEL=false \
    POETRY_CACHE_DIR=/poetrycache/ \
    NUMBA_CACHE_DIR=/tmp

RUN apt-get update && apt-get install -y \
    python3.8 \
    python3.8-dev \
    python3-pip \
    libxml2-dev \
    libgd-perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libproj-dev \
 && pip3 install poetry \
 && ln -s /usr/bin/python3.8 /usr/bin/python

RUN install2.r --error --deps TRUE -r $CRAN --skipinstalled \
    argparse \
    cowplot \
    digest \
    dplyr \
    dynamicTreeCut \
    flashClust \
    forcats \
    ggnewscale \
    ggplot2 \
    ggraph \
    ggrepel \
    ggupset \
    gprofiler2 \
    gtools \
    here \
    jsonlite \
    kableExtra \
    matrixStats \
    missForest \
    naniar \
    pheatmap \
    purrr \
    RColorBrewer \
    readr \
    remotes \
    tibble \
    tidyr \
    UpSetR \
    yaml \
    littler \
    tidyverse \
    msigdbr \
  && /usr/local/lib/R/site-library/littler/examples/installBioc.r \
    biomaRt \
    limma \
    DESeq2 \
    STRINGdb \
    org.Mm.eg.db \
    org.Hs.eg.db \
    impute \
    GO.db \
    AnnotationDbi \
    EnhancedVolcano \
    GeneTonic \
    GSVA \
    proDA \
  && Rscript -e 'BiocManager::install("preprocessCore", configure.args="--disable-threading")' \
  && rm -rf /tmp/downloaded_packages

#protti needs STRINGdb, limma and proDA installed first
RUN install2.r --error --deps TRUE -r $CRAN --skipinstalled \
    protti

RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r \
  WGCNA \
  clusterProfiler \
  enrichplot \
  rtracklayer \
  pathview \
  DOSE \
  DRIMSeq\
  DEXSeq \
  stageR\
  tximport\
  GenomicFeatures\
  topGO \
  pcaExplorer \
  && rm -rf /tmp/downloaded_packages

#https://github.com/hms-dbmi/UpSetR/pull/214
RUN Rscript -e 'remotes::install_github("zeehio/UpSetR@fix-single-element")'

RUN Rscript -e 'remotes::install_github("ctlab/fgsea@v1.19.2")'

# Let rpy2 find local R installation
RUN echo "/usr/local/lib/R/lib/" >> /etc/ld.so.conf.d/local-R-install.conf && ldconfig

RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 \
    && mkdir -p samtools \
    && tar -xf samtools-1.14.tar.bz2 -C samtools --strip-components=1 \
    && cd samtools \
    && ./configure --without-curses --prefix=/usr/local/ \
    && make install \
    && make -C htslib-1.14 install \
    && cd .. \
    && rm -r samtools samtools-1.14.tar.bz2
RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 \
    && mkdir -p bcftools \
    && tar -xf bcftools-1.14.tar.bz2 -C bcftools --strip-components=1 \
    && cd bcftools \
    && ./configure --prefix=/usr/local/ \
    && make install \
    && cd .. \
    && rm -r bcftools bcftools-1.14.tar.bz2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip \
    && unzip bowtie2-2.4.4-linux-x86_64.zip \
    && cp bowtie2-2.4.4-linux-x86_64/bowtie* /usr/local/bin \
    && rm -r bowtie2-2.4.4-linux-x86_64/ bowtie2-2.4.4-linux-x86_64.zip \
    && wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.0.tar.gz \
    && tar -xzf v0.15.0.tar.gz \
    && cp -r FastQ-Screen-0.15.0/* /usr/local/bin \
    && rm -r FastQ-Screen-0.15.0/ v0.15.0.tar.gz \
    && perl -MCPAN -e 'install GD::Graph'

COPY pyproject.toml poetry.lock ./

RUN mkdir -p /poetrycache/ /.config/ /.cache/ && chmod uog+rw -R /poetrycache/ /.config/ /.cache/ && poetry install --no-root

ENTRYPOINT ["/usr/local/bin/poetry", "run"]
