# Use the official Ubuntu 24.04 as the base image
FROM ubuntu:24.04

# Arguments for non-root user
ARG USERNAME=devuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Europe/Paris

# Install dependencies and tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    ca-certificates \
    sudo \
    python3 \
    python3-pip \
    python3-venv \
    libgtest-dev \
    libeigen3-dev \
    libboost-dev \
    && rm -rf /var/lib/apt/lists/*

# Setup Python virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install Python packages and gitman
RUN pip install --no-cache-dir \
    pre-commit \
    cpplint \
    black \
    flake8 \
    gitman \
    numpy \
    matplotlib \
    scipy \
    && ln -s /opt/venv/bin/gitman /usr/local/bin/gitman \
    && ln -s /opt/venv/bin/pre-commit /usr/local/bin/pre-commit

# Copy Gitman config and install dependencies



# Setup non-root user
# Ubuntu 24.04 comes with a default 'ubuntu' user (UID 1000). We rename it to USERNAME.
RUN usermod -l $USERNAME ubuntu \
    && usermod -d /home/$USERNAME -m $USERNAME \
    && groupmod -n $USERNAME ubuntu \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set working directory and ownership
WORKDIR /app/MPCC
RUN chown -R $USERNAME:$USERNAME /app && \
    chown -R $USERNAME:$USERNAME $VIRTUAL_ENV

# Switch to non-root user
# USER $USERNAME

# Configure Git to avoid "dubious ownership" errors for the mounted volume
RUN git config --global --add safe.directory /app/MPCC

# Define volume
VOLUME /app/MPCC

CMD ["/bin/bash"]
