# Base Image
FROM assemblyrepair:latest

COPY src/ /workdir/
RUN chmod +x /workdir/*.sh

WORKDIR /workdir/

CMD ["/bin/bash"]
