FROM haoyangz/keras-docker
MAINTAINER Haoyang Zeng  <haoyangz@mit.edu>

RUN pip install --upgrade --no-deps git+https://github.com/maxpumperla/hyperas@0.1.2
RUN pip install --upgrade --no-deps hyperopt pymongo scikit-learn networkx

ENV THEANO_FLAGS='cuda.root=/usr/local/cuda,device=gpu0,floatX=float32,lib.cnmem=0.8'

# Install CUDA repo (needed for cuDNN)
ENV CUDA_REPO_PKG=cuda-repo-ubuntu1404_7.0-28_amd64.deb
RUN wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1404/x86_64/$CUDA_REPO_PKG && \
  dpkg -i $CUDA_REPO_PKG

# Install cuDNN v4
ENV ML_REPO_PKG=nvidia-machine-learning-repo_4.0-2_amd64.deb
RUN wget http://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu1404/x86_64/$ML_REPO_PKG && \
  dpkg -i $ML_REPO_PKG && \
    apt-get update && apt-get install -y libcudnn4 libcudnn4-dev

COPY main.py /scripts/
COPY cnn /scripts/cnn/
COPY helper /scripts/helper/
RUN cd /scripts/;wget http://gerv.csail.mit.edu/MethylDecoder_models.tar.gz -q;tar -zxvf MethylDecoder_models.tar.gz
WORKDIR /scripts/
