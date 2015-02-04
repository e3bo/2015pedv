FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --error \
  fields \
  ggplot2 \
  igraph \
  knitr \
  pander \
  vegan \
  && rm -rf /tmp/download_packages/ /tmp/*.rds

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  poppler-utils \
  python3-lxml \
  python3-pkg-resources \
  python3-pandas

CMD ["bash"]
