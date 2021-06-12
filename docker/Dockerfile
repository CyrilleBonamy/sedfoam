FROM openfoamplus/of_v2006_centos73
LABEL maintainer "cyrille.bonamy@univ-grenoble-alpes.fr"
ARG WM_NCOMPPROCS=10

RUN yum update -y \
  && yum install -y \
  git python3 python3-pip unzip \
  mercurial readline-devel \
  && yum clean all

RUN pip3 install fluidfoam


USER 1000

WORKDIR /home/ofuser
RUN /bin/bash -c "hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam -r v2020.06"

WORKDIR /home/ofuser/swak4foam
RUN /bin/bash -c "shopt -s expand_aliases && source /opt/OpenFOAM/setImage_v2006.sh && ./maintainanceScripts/compileRequirements.sh"
RUN /bin/bash -c "shopt -s expand_aliases && source /opt/OpenFOAM/setImage_v2006.sh && ./Allwmake"

WORKDIR /home/ofuser
RUN rm -rf swak4foam
#RUN git clone --recurse-submodules -b module_openfoam.com https://github.com/SedFoam/sedfoam.git
RUN git clone --recurse-submodules https://github.com/SedFoam/sedfoam.git

WORKDIR /home/ofuser/sedfoam
RUN /bin/bash -c "shopt -s expand_aliases && source /opt/OpenFOAM/setImage_v2006.sh && ./Allwmake && ./Allclean"
#RUN find /home/ofuser/sedfoam -name "*.o" -exec rm -f {} \;

#WORKDIR /opt
#WORKDIR /opt/swak4Foam
#RUN /bin/bash -c "source /opt/openfoam/5.0/OpenFOAM-5.0/etc/bashrc && ./maintainanceScripts/compileRequirements.sh"
##COPY /opt/swak4Foam/swakConfiguration.debian  /opt/swak4Foam/swakConfiguration
#RUN /bin/bash -c "source /opt/openfoam/5.0/OpenFOAM-5.0/etc/bashrc && ./Allwmake"
#RUN /bin/bash -c "source /opt/openfoam/5.0/OpenFOAM-5.0/etc/bashrc && ./maintainanceScripts/copySwakFilesToSite.sh"

WORKDIR /home/ofuser
ENV HOME /home/ofuser
ENTRYPOINT ["/bin/bash", "-c"]
