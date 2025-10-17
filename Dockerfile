FROM pytorch/pytorch:2.8.0-cuda12.6-cudnn9-runtime

# Update and install dependencies
USER root
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install gcc g++ git -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy DECIPHER package and install
WORKDIR /app
COPY ./decipher ./decipher
COPY ./README.md ./README.md
COPY ./pyproject.toml ./pyproject.toml
RUN pip --no-cache-dir install -e "." && rm -rf /tmp/* && pip cache purge && install_pyg_dependencies && rm -rf /tmp/* && pip cache purge
