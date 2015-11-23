FROM r-base
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  ghostscript \
  libgdal1-dev \
  libgfortran-5-dev \ 
# ^needed for some R packages to build
  libproj-dev \
  openssh-server \
  poppler-utils \
  python3-lxml \
  python3-pkg-resources \
  python3-pandas \
  r-cran-car \
  r-cran-coda \
  r-cran-doparallel \
  r-cran-hmisc \
  r-cran-ggplot2 \
  r-cran-igraph \
  r-cran-lme4 \
  r-cran-maps \
  r-cran-mapproj \
  r-cran-maptools \
  r-cran-polspline \
  r-cran-plyr \
  r-cran-pscl \
  r-cran-rcolorbrewer \
  r-cran-reshape \
  r-cran-reshape2 \
  r-cran-rgl \
  r-cran-scales \
  r-cran-sp \
  r-cran-surveillance \
  r-cran-vcd \
  r-cran-vegan \
  r-cran-xml

RUN install2.r --error \
  c060 \
  DiceDesign \ 
  DiceEval \ 
  DiceKriging \
  DiceView \
  DescTools \
  fields \
  geoR \
  GGally \
  glmnet \
  gridBase \
  grImport \
  knitr \
  lhs \
  mda \
  pander \
  randtoolbox \
  raster \
  rgdal \
  sensitivity \
  && rm -rf /tmp/download_packages/ /tmp/*.rds

RUN install2.r --error R2admb \
  && install2.r --repos http://glmmadmb.r-forge.r-project.org/repos --error glmmADMB \
  && install2.r --repos http://www.math.mcmaster.ca/bolker/R --error coefplot2 \
  && rm -rf /tmp/download_packages/ /tmp/*.rds

RUN mkdir /var/run/sshd && echo 'docker:screencast' | chpasswd \
  && sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
EXPOSE 22

CMD ["/usr/sbin/sshd", "-D"]
