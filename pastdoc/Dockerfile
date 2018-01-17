FROM ubuntu:xenial

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# Create the development environment
RUN apt update && \
    apt install -y python3-setuptools python3-dev python3-pip python3
RUN pip3 install --no-cache-dir --upgrade pip

# Copy in the main code
ADD cytopast_main.py /
# Install cytopast
RUN pip3 install --no-cache-dir cytopast

# Clean up the development environment to reduce image size
RUN apt remove -y build-essential python3.*-dev && \
    apt autoremove -y

# Copy in PASTML and make it executable
ADD https://github.com/saishikawa/PASTML/releases/download/0.1/PASTML /apps/
RUN chmod +x /apps/PASTML
ENV PATH="/apps/:${PATH}"

WORKDIR /tmp

# The entrypoint runs cytopast with command line arguments
ENTRYPOINT ["/usr/bin/python3", "/cytopast_main.py"]
