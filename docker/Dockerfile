#
# https://github.com/cpp-projects-showcase/docker-images/tree/master/ubuntu2004
#
FROM docker.io/ubuntu:20.04

# Environment
ENV container docker
ENV HOME /home/build
ENV LANGUAGE en_US:en
ENV LANG en_US.UTF-8
ENV LC_ALL $LANG

# Update the system
RUN echo "Europe/Berlin" > /etc/timezone
RUN echo 'APT::Get::Assume-Yes "true";' > /etc/apt/apt.conf.d/90-yes
RUN apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends \
	apt-utils apt-transport-https \
	ca-certificates locales tzdata sudo \
	zip unzip gzip bzip2 xz-utils tar unrar p7zip-full \
	curl wget netcat net-tools aptitude

# Generate the locales
RUN locale-gen $LANG && update-locale LANG=$LANG LANGUAGE=$LANGUAGE LC_ALL=$LC_ALL

# Configure the time-zone
RUN dpkg-reconfigure -f noninteractive tzdata

# Basic, C++ and Python packages
RUN apt-get -y install less htop screen \
	git keychain gawk \
	bash-completion vim-nox emacs-nox apt-utils keyutils ftp \
	zlib1g-dev libbz2-dev \
	lsb-release libgmp-dev \
	gcc g++ cppcheck clang cmake manpages patch pkg-config \
	m4 autoconf automake libtool libltdl-dev build-essential \
	flex bison \
	libboost-all-dev libxapian-dev \
	libreadline-dev libncurses5-dev \
	libzmqpp-dev libssl-dev libffi-dev \
	swig graphviz libopenblas-dev && \
    apt-get -y install libmpich-dev libopenmpi-dev \
	sqlite3 libsqlite3-dev \
	mysql-client libmysqlclient-dev libmysql++-dev \
	postgresql-server-dev-all \
	libicu-dev libprotobuf-dev protobuf-compiler \
	python2 libpython2-dev python3 python3-dev libpython3-dev python3-pip \
	python3-django libapache2-mod-wsgi-py3 \
	libgeos++-dev \
	doxygen ghostscript texlive-latex-recommended \
	r-base r-base-dev \
	rake

RUN apt-get update \
	&& apt-get -y install build-essential \
	&& apt-get install -y wget \
	&& rm -rf /var/lib/apt/lists/* \
	&& wget https://github.com/Kitware/CMake/releases/download/v3.24.1/cmake-3.24.1-Linux-x86_64.sh \
		-q -O /tmp/cmake-install.sh \
		&& chmod u+x /tmp/cmake-install.sh \
		&& mkdir /opt/cmake-3.24.1 \
		&& /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake-3.24.1 \
		&& rm /tmp/cmake-install.sh \
		&& ln -s /opt/cmake-3.24.1/bin/* /usr/local/bin

# Install the packages for Python development

# RUN apt-get update && \
# 	apt-get install -y python3 python3-dev python3-pip && \

RUN pip3 install --no-cache-dir "pybind11[global]" && \
	pip3 install --no-cache-dir ipykernel && \
	pip3 install --no-cache-dir numpy matplotlib 

# Cleaning
RUN apt-get clean && rm -rf /tmp/* /var/tmp/*

# Create the `build` user (for the development activities)
RUN adduser --disabled-password --gecos "" build
RUN echo "build ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/build && \
    chmod 0440 /etc/sudoers.d/build

# Configure SSH
RUN mkdir -p $HOME/.ssh && chmod 700 $HOME/.ssh

# Switch to the 'build' user
WORKDIR $HOME
USER build

# Git prompt
#RUN git clone https://github.com/magicmonty/bash-git-prompt.git ~/.bash-git-prompt --depth=1


# RUN git clone https://github.com/pyenv/pyenv.git $HOME/.pyenv
# ENV PATH $HOME/.pyenv/bin:$HOME/.pyenv/shims:$PATH
# RUN pyenv install 3.9-dev && \
#     pyenv global 3.9 && \
#     pip install -U pip pipenv && \
# 	pip install pybind11
# RUN pyenv global system || echo "No default system version of Python. Sticking to 3.9.0"


# Entry point
CMD ["/bin/bash"]
