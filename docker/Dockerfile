# Use identical base image (20.04) to what we use in the fv3net prognostic run:
# https://github.com/ai2cm/fv3net/blob/382409f76fd26d70457b349a48fa0e63f4f26303/Makefile#L12
#
# This is convenient for a couple reasons:
#
# * It installs Python 3.8 as the system Python, which should allow us to
#   seamlessly install the additional Python packages used in the existing
#   prognostic_run image, to their exact versions.
#
# * It installs a version of MPICH that does not enable link time optimization
#   (LTO) by default.  Newer Ubuntu LTS base images, e.g. 22.04, install
#   versions of MPICH that enable link time optimization by default in their
#   compiler wrappers (i.e. mpif90 and the like; see the presence of -flto=auto
#   and other LTO-related flags in the output of `mpif90 -show`).  LTO can lead
#   to differences in answers between applications that were linked differently
#   (e.g. the fortran executable and the Python-wrapped model). Compiling
#   without LTO allows the fortran executable and the Python wrapper to
#   robustly produce bitwise identical results.  See
#   https://github.com/ai2cm/SHiELD-minimal/pull/6 for some details on MPICH in
#   the Ubuntu LTS 22.04 image.  This can be addressed by manually adding
#   `-fno-lto` compiler flags to Makefile templates, but it is easiest not to 
#   need to worry about this.
FROM ubuntu@sha256:9101220a875cee98b016668342c489ff0674f247f6ca20dfc91b91c0f28581ae

# Install and use GNU version 10 compiliers instead of the default version 9,
# since the latest version of FMS fails to build with GNU version 9
# compilers.
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install \
    cmake \
    gcc-10 \
    gfortran-10 \
    g++-10 \
    git \
    libmpich-dev \
    libnetcdf-dev \
    libnetcdff-dev \
    m4 \
    make \
    mpich \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-setuptools

# GNU version 9 compilers are unavoidably installed by default; using
# update-alternatives ensures that gcc, gfortran, and g++ point to GNU version
# 10 compilers, including within MPI-wrapped compilers.
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 1
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-10 1
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 1

ENV SHiELD_FC=mpif90
ENV SHiELD_CC=mpicc
ENV SHiELD_CXX=mpicxx
ENV SHiELD_LD=mpif90

ENV SUBMODULE_DIR=/SHiELD
COPY submodules ${SUBMODULE_DIR}/

# Build FMS, NCEPlibs, and SHiELD using SHiELD_build with position independent
# code FMS_CPPDEFS is needed to address
# https://github.com/NOAA-GFDL/FMS/issues/426
RUN cd ${SUBMODULE_DIR}/SHiELD_build/Build && \
    FMS_CPPDEFS="-DHAVE_GETTID" \
    FC=${SHiELD_FC} \
    CC=${SHiELD_CC} \
    CXX=${SHiELD_CXX} \
    LD=${SHiELD_LD} \
    TEMPLATE=site/gnu.mk \
    AVX_LEVEL=-march=native \
    ./COMPILE shield 64bit gnu pic

COPY requirements.txt /tmp
RUN \
    FC=${SHiELD_FC} \
    CC=${SHiELD_CC} \
    CXX=${SHiELD_CXX} \
    LD=${SHiELD_LD} \
    python3 -m pip install -r /tmp/requirements.txt

# Prescribe LDSHARED as well to ensure that the wrapper is built with the GNU
# 10 compilers.  See https://github.com/ai2cm/fv3gfs-fortran/issues/330 for
# more comprehensive discussion.  This may not be absolutely necessary (we do
# not worry about this in fv3net), but it feels cleaner not to mix compilers if
# we can avoid it.
COPY wrapper /wrapper
RUN make -C wrapper \
    FC=${SHiELD_FC} \
    CC=${SHiELD_CC} \
    CXX=${SHiELD_CXX} \
    LD=${SHiELD_LD} \
    LDSHARED="${SHiELD_CC} -shared" \
    PIC=Y \
    OPENMP=Y \
    AVX=Y \
    build install
RUN python3 -c "import shield.wrapper"
