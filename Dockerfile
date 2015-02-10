FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  poppler-utils \
  python3-lxml \
  python3-pkg-resources \
  python3-pandas

RUN install2.r --error \
  car \
  coda \
  c060 \
  fields \
  geoR \
  ggplot2 \
  glmnet \
  Hmisc \
  igraph \
  knitr \
  lme4 \
  maps \
  pander \
  plyr \
  pscl \
  reshape \
  vegan \
&& rm -rf /tmp/download_packages/ /tmp/*.rds
RUN install2.r --repos http://r-forge.r-project.org --error glmmADMB \
&& install2.r --repos http://www.math.mcmaster.ca/bolker/R --error coefplot2 \
&& rm -rf /tmp/download_packages/ /tmp/*.rds

CMD ["bash"]
