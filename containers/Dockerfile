ARG ARG_PYTHON_VERSION=3.11
ARG ARG_RUNTIME_IMAGE=${ARG_PYTHON_VERSION}-bookworm

FROM python:${ARG_RUNTIME_IMAGE}

# install dependencies
RUN apt-get update && apt-get install -y libblas-dev liblapack-dev gfortran r-base
RUN apt-get autoremove && apt-get clean
RUN apt-get install -y --no-install-recommends git 
RUN pip install mkl

# set working directory
WORKDIR /home

# copy over python requirements
COPY requirements.txt requirements.txt

# install all python requirements
RUN pip install --no-cache-dir -r requirements.txt

# install required R packages
RUN R -e "install.packages('sctransform', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install_version('sctransform', version='0.4.1', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# set entry point as a bash command
CMD ["bash"]
