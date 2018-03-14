#!/bin/bash

DATABASE=
QUERY=
FORMAT="docsum"
FILE_STRING="_genomic.fna.gz"

while [ "$1" != "" ]; do
    case $1 in
        -db | --database )      shift
                                DATABASE=$1
                                ;;
        -q | --query )          shift
                                QUERY=$1
                                ;;
        -f | --format )          shift
                                FORMAT=$1
                                ;;
        -fs | --file-string )   shift
                                FILE_STRING=$1
                                ;;
        * )                     echo "Command not found"
                                exit 1
    esac
    shift
done

wget `esearch -db "$DATABASE" -query "$QUERY" | \
efetch -format "$FORMAT" | \
xtract -pattern DocumentSummary -element FtpPath_GenBank | \
awk -F"/" '{print $0"/"$NF"'$FILE_STRING'"}'`
