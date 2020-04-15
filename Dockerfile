FROM nfcore/base:1.9
LABEL authors="Francesco Lescai" \
      description="Docker image containing all software requirements for the nibscbioinformatics/scoop pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nibscbioinformatics-scoop-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nibscbioinformatics-scoop-1.0dev > nibscbioinformatics-scoop-1.0dev.yml

# Both Humann2 and Metaphlan2 are a difficult with databases
# So packing default databases within the container to simplify execution
# Although this will make container larger in size

RUN mkdir -p /opt/databases/humann2_dbs
RUN cd /opt/databases
RUN humann2_databases --download chocophlan full humann2_dbs
RUN humann2_databases --download uniref uniref50_ec_filtered_diamond humann2_dbs
RUN humann2_databases --download uniref uniref50_GO_filtered_rapsearch2 humann2_dbs
RUN humann2_databases --download uniref uniref50_diamond humann2_dbs
RUN humann2_databases --download uniref uniref90_ec_filtered_diamond humann2_dbs
RUN humann2_databases --download utility_mapping full humann2_dbs

# Metaphlan has its own way to handle this so we have to run a short analysis to trigger download

RUN mkdir -p /opt/tests
RUN cd /opt/tests
RUN wget https://bitbucket.org/biobakery/metaphlan2/downloads/SRS019033.fastq
RUN metaphlan2.py SRS019033.fastq --input_type fastq > profiled_metagenome.txt
