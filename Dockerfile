FROM nfcore/base:1.12.1
LABEL authors="Sabrina Krakau and Leon Kuchenbecker" \
    description="Docker image containing all software requirements for the nf-core/metapep pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a
# Manually install epytope / FRED-2 from GH
RUN /opt/conda/envs/nf-core-metapep-1.0dev/bin/pip install git+https://github.com/KohlbacherLab/epytope.git@a863afc5131a33d3510ba0a397cd34bbcaae6270

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-metapep-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-metapep-1.0dev > nf-core-metapep-1.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

# Download MHCflurry models
ENV MHCFLURRY_DATA_DIR /mhcflurry-data
ENV MHCFLURRY_DOWNLOADS_CURRENT_RELEASE 1.4.0
RUN mkdir -p "$MHCFLURRY_DATA_DIR"
RUN mhcflurry-downloads fetch models_class1
