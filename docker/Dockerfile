FROM ubuntu:xenial
LABEL maintainer="ekkus93@gmail.com"

RUN apt-get update && apt-get upgrade -y 
RUN apt-get install wget screen bzip2 zlib1g-dev unzip curl -y

RUN apt-get install cmake make build-essential git -y

RUN apt-get install openssl libssl-dev libuv1.dev -y

# common
RUN rm -f /usr/lib/libuWS.so
RUN cd /tmp && git clone https://github.com/udacity/CarND-MPC-Project.git
COPY ./install-ubuntu.sh /tmp/CarND-MPC-Project/install-ubuntu.sh
RUN cd /tmp/CarND-MPC-Project && sh ./install-ubuntu.sh

# ipopt
RUN apt-get install gfortran -y
RUN cd /tmp && wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.7.zip && unzip Ipopt-3.12.7.zip && rm Ipopt-3.12.7.zip
COPY ./install_ipopt.sh /tmp/Ipopt-3.12.7
RUN cd /tmp/Ipopt-3.12.7 && bash ./install_ipopt.sh Ipopt-3.12.7

# CppAD
RUN apt-get install cppad -y

# create ubuntu user
RUN useradd -ms /bin/bash ubuntu
RUN mkdir -p /home/ubuntu/website && chown -R ubuntu /home/ubuntu

# set bash shell
RUN chsh -s /bin/bash ubuntu

# set home dir for ubuntu
RUN usermod -d /home/ubuntu ubuntu

VOLUME ["/home/ubuntu/work"]

COPY ./version.txt /tmp/

USER ubuntu
WORKDIR "/home/ubuntu/work"

CMD ["bash"]

EXPOSE 4567

