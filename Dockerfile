FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --error \
  car \
  coda \
  fields \
  ggplot2 \
  igraph \
  knitr \
  lme4 \
  pander \
  pscl \
  vegan \
  && rm -rf /tmp/download_packages/ /tmp/*.rds
RUN install2.r --repos http://r-forge.r-project.org --error glmmADMB \
  && rm -rf /tmp/download_packages/ /tmp/*.rds

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  poppler-utils \
  python3-lxml \
  python3-pkg-resources \
  python3-pandas

CMD ["bash"]
