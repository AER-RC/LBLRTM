# Build LBLRTM and LBLRTM in the latest CentOS environment
FROM centos:6

# install tools necessary for building models
RUN yum -y install make gcc-gfortran && \
  yum -y -q update && \
  yum clean all
RUN rm -rf /var/cache/yum/*

WORKDIR /LBLRTM

# necessary code for LBLRTM
ADD aer_rt_utils /LBLRTM/aer_rt_utils
ADD build /LBLRTM/build
ADD src /LBLRTM/src

# build the model and clean up after the build
RUN cd /LBLRTM/build; \
  make -f make_lblrtm linuxGNUdbl; \
  rm -rf lblrtm_v12.9_linux_gnu_dbl.obj *.mod; cd /LBLRTM; \
  ln -s lblrtm_v12.9_linux_gnu_dbl lblrtm

VOLUME LBLRTM_In
VOLUME LBLRTM_Out

COPY LBLRTM_entrypoint.sh .

# run model and build binary line file (TAPE3)
ENTRYPOINT ["./LBLRTM_entrypoint.sh"]
