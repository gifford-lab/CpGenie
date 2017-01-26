FROM haoyangz/keras-docker:cuda7.0
MAINTAINER Haoyang Zeng  <haoyangz@mit.edu>

RUN pip install --upgrade --no-deps git+https://github.com/maxpumperla/hyperas@0.1.2
RUN pip install --upgrade --no-deps hyperopt pymongo scikit-learn networkx

ENV THEANO_FLAGS='cuda.root=/usr/local/cuda,device=gpu0,floatX=float32,lib.cnmem=0.1,base_compiledir=/runtheano/.theano'

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('snow')"
RUN apt-get update;apt-get -y install lzop

COPY main.py /scripts/
COPY variant.py /scripts/
COPY cnn /scripts/cnn/
COPY helper /scripts/helper/
COPY data /scripts/data/
RUN cd /scripts/;wget http://gerv.csail.mit.edu/CpGenie_models.tar.gz -q;tar -zxvf CpGenie_models.tar.gz
RUN cd /scripts/data;wget http://gerv.csail.mit.edu/hg19.in.lzo -q; lzop -d hg19.in.lzo
RUN mkdir /runtheano/;chmod -R 777 /runtheano/
WORKDIR /scripts/
