# Use Ubuntu 22.04 which includes Python 3.10
FROM ubuntu:22.04

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Update apt and install required packages:
# - cmake, clang, llvm, ninja-build, git for building
# - python3.10 and related packages for Python
# - lcov for coverage tools
RUN apt update && apt install -y \
    git \
    clang \
    cmake \
    libopenblas-dev \
    liblapack-dev \
    lcov \
    llvm \
    python3.10 \
    python3.10-dev \
    python3.10-venv && \
    rm -rf /var/lib/apt/lists/*

# Set PATH to include /usr/local/bin
ENV PATH="/usr/local/bin:${PATH}"

# Copy the requirements.txt file
COPY requirements.txt .

# Create virtual environment and install dependencies
RUN python3 -m venv /venv && \
    /venv/bin/pip install --upgrade pip && \
    /venv/bin/pip install -r requirements.txt

# Add the virtual environment's bin directory to PATH.
ENV PATH="/venv/bin:${PATH}"

# By default, run a shell.
CMD ["/bin/bash"]
