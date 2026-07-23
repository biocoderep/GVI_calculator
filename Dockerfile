FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    curl \
    default-jre \
    libhmsbeagle-dev \
    libhmsbeagle-java \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /app

# Download and install IQ-TREE
RUN wget -qO iqtree.tar.gz https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz && \
    tar -xzf iqtree.tar.gz && \
    cp iqtree-2.3.6-Linux-intel/bin/iqtree2 /usr/local/bin/iqtree2 && \
    rm -rf iqtree.tar.gz iqtree-2.3.6-Linux-intel

# Install BEAST2 and missing packages
RUN wget -qO beast.tgz https://github.com/CompEvol/beast2/releases/download/v2.6.7/BEAST.v2.6.7.Linux.tgz && \
    tar -xzf beast.tgz && \
    mv beast /usr/local/beast && \
    ln -s /usr/local/beast/bin/beast /usr/local/bin/beast && \
    rm beast.tgz

# Copy requirements and install
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application
COPY . .

# Set Environment variables
ENV DOCKER_ENV=true
ENV FLASK_APP=app.py
ENV FLASK_RUN_HOST=0.0.0.0

# Expose port
EXPOSE 5000

# Start the Flask app
CMD ["flask", "run"]
