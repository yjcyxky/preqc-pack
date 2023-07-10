FROM rust:1.54.0-buster
LABEL Build environment for musl-demo project

RUN echo 'deb http://ftp.cn.debian.org/debian/ stable main' > /etc/apt/sources.list
RUN apt update
RUN apt install libssl-dev libsodium-dev build-essential binutils upx -y

CMD ["/bin/bash"]