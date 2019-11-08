FROM rocker/shiny

LABEL maintainer "Keqiang Li <keqianglli@gmail.com>"

RUN R -e "install.packages(c('shinythemes', 'DT', 'remotes'), repos='https://cloud.r-project.org/')"

RUN R -e "install.packages(c('magrittr', 'purrr', 'rlang', 'glue', 'dplyr', 'readr'), repos='https://cloud.r-project.org/')"

RUN R -e "remotes::install_github(c('keqiang/shinyfio', 'keqiang/genemap'))"

WORKDIR /srv/shiny-server/genemap

ADD *.R ./
ADD sample_data ./sample_data
