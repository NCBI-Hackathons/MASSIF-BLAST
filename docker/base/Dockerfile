# Base Image
FROM centos:centos7

RUN yum install -y \
    wget \
    perl-HTML-Parser \
    # cpan \
    perl-libwww-perl \
    # perl-Crypt-SSLeay \
    perl-LWP-Protocol-https \
    java-1.8.0-openjdk-devel \
    zlib-devel

# Install sbt
RUN curl https://bintray.com/sbt/rpm/rpm | tee /etc/yum.repos.d/bintray-sbt-rpm.repo && \
    yum install -y sbt

# Needed for edirect
RUN perl -MCPAN -e 'install Bundle::LWP'

# Install edirect
RUN cd /lib/ && \
    /bin/bash && \
    perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); \
        $ftp->login; $ftp->binary; \
        $ftp->get("/entrez/entrezdirect/edirect.tar.gz");' \
    gunzip -c edirect.tar.gz | tar xf - && \
    rm edirect.tar.gz && \
    builtin exit && \
    export PATH=${PATH}:/lib/edirect >& /dev/null || setenv PATH "${PATH}:/lib/edirect" && \
    ./edirect/setup.sh

# Install Magic-BLAST
RUN cd /lib/ && \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/1.3.0/ncbi-magicblast-1.3.0-x64-linux.tar.gz && \
    tar -xzf ncbi-magicblast-1.3.0-x64-linux.tar.gz && \
    rm ncbi-magicblast-1.3.0-x64-linux.tar.gz

# Install BLAST
RUN cd /lib/ && \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz

# Install Pilon
RUN cd /lib/ && \
    wget https://github.com/broadinstitute/pilon/releases/download/v1.22/pilon-1.22.jar

# Install samtools
RUN cd /lib/ && \
    wget -O samtools.tar.gz https://github.com/samtools/samtools/archive/1.7.tar.gz && \
    tar -xzf samtools.tar.gz && \
    rm samtools.tar.gz && \
    cd samtools/ && \
    ./configure && \
    make && \
    make install

# Get FALP executable
ADD falp_1.06/falp.exe.gz /lib/falp.exe.gz
RUN gunzip /lib/FALP_1.06.zip

USER dockeruser
WORKDIR /home/dockeruser

# Add edirect to the Path
ENV PATH "$PATH:/lib/edirect"
# Add Magic-BLAST to the Path
ENV PATH "$PATH:/lib/ncbi-magicblast-1.3.0/bin"
# Add BLAST to the Path
ENV PATH "$PATH:/lib/ncbi-blast-2.7.1+/bin"

# Create Alias to run Pilon
RUN echo 'alias pilon="/lib/pilon-1.22.jar"' >> ~/.bashrc
# Create Alias to run FALP
RUN echo 'alias falp="/lib/falp.exe"' >> ~/.bashrc


CMD ["/bin/bash"]
