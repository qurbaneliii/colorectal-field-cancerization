FROM rocker/r-ver:4.5.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv make git libcurl4-openssl-dev \
    libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libgit2-dev \
    libglpk-dev libgsl-dev pandoc \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /analysis

# Restore dependency graphs before copying frequently changing source files.
COPY renv.lock pyproject.toml requirements-lock.txt ./
RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org')" \
    && Rscript -e "renv::restore(prompt = FALSE)"
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --requirement requirements-lock.txt \
    && /opt/venv/bin/pip check
ENV PATH="/opt/venv/bin:${PATH}"

COPY . .
RUN /opt/venv/bin/pip install --no-cache-dir --no-deps -e . \
    && Rscript R/verify_environment.R \
    && python -c "import numpy, pandas, sklearn, pyarrow, yaml"

# The image verifies environments only. Mount data and choose an explicit runtime target.
CMD ["make", "help"]
