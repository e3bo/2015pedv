FROM r-base:3.1.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  ghostscript \
  poppler-utils \
  python3-lxml \
  python3-pkg-resources \
  python3-pandas \
  r-cran-xml

RUN install2.r --error \
  BatchExperiments \
  car \
  coda \
  c060 \
  DescTools \
  fields \
  geoR \
  GGally \
  ggplot2 \
  glmnet \
  gridBase \
  grImport \
  Hmisc \
  igraph \
  knitr \
  lme4 \
  mapproj \
  maps \
  pander \
  plyr \
  pscl \
  reshape \
  reshape2 \
  RColorBrewer \
  scales \
  vcd \
  vegan \
&& rm -rf /tmp/download_packages/ /tmp/*.rds
RUN install2.r --repos http://r-forge.r-project.org --error glmmADMB \
&& install2.r --repos http://www.math.mcmaster.ca/bolker/R --error coefplot2 \
&& rm -rf /tmp/download_packages/ /tmp/*.rds

CMD ["bash"]
