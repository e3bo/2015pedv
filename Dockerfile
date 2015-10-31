FROM r-base
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  ghostscript \
  libgdal1 \
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
  maptools \
  pander \
  plyr \
  pscl \
  RColorBrewer \
  reshape \
  reshape2 \
  rgdal \
  scales \
  surveillance \
  vcd \
  vegan \
&& rm -rf /tmp/download_packages/ /tmp/*.rds

RUN install2.r --error R2admb
&& install2.r --repos http://glmmadmb.r-forge.r-project.org/repos --error glmmADMB \
&& install2.r --repos http://www.math.mcmaster.ca/bolker/R --error coefplot2 \
&& rm -rf /tmp/download_packages/ /tmp/*.rds

CMD ["bash"]
