# Base Docker image for project

## Build base image

```
./build.sh
```

## FALP

The intent is not to redistribute others software, however the zip file provided (here)[ftp://ftp.ncbi.nih.gov/pub/spouge/web/software/FALP_1.06/FALP_1.06.zip] requires zip version >= 2.0. The base image, CentOS 7, has `gzip` version 1.5.

## Running software in Docker

Commands to run tools installed in Docker image.

### EDirect Command Line tools

```
docker run assemblyrepair:latest esearch --help
```

### Magic-BLAST

```
docker run assemblyrepair:latest magicblast --help
```

### BLAST

```
docker run assemblyrepair:latest blastn --help
```

### Pilon

```
docker run assemblyrepair:latest java -jar /lib/pilon-1.22.jar --help
```

### Samtools

```
docker run assemblyrepair:latest samtools --help
```

### FALP

```
docker run assemblyrepair:latest /lib/falp.exe -help
```
