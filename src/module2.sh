#!/bin/bash

## MODULE to run Pilon to fix assembly using DNA data
## written by Kyle Roell

#set up directory structure for code
DIR=$(mktemp -d)
trap "rm -rf $DIR" EXIT

acc=false
dnafile=false
keep=false

size=5

#rm input_dna.out 2>/dev/null
#rm input_rna.out 2>/dev/null

#########################
# The command line help1 #
#########################
display_help() {
printf "\nUSAGE\n\n"
printf "   -h, --help                 Display this dialog\n\n"
printf "   --genome                   Specify genome assembly data to be used\n"
printf "                              This is MANDATORY!\n\n"
printf "   --DNAfiles                 Specify SRA DNA fastq input files to be used with Pilon,\n"
printf "                              multiple files may be input separated by spaces\n\n"
printf "   --acc, --accession         Specify SRA DNA accession numbers to be used with Pilon,\n"
printf "                              multiple numbers may be input separated by spaces\n\n"
printf "   --organism                 Specify organism to be used in NCBI SRA data search,\n"
printf "                              must be the same as organism used in --genome\n"
printf "                              This only applies when no SRAs are input using --accesion or --DNAfiles\n\n"
printf "   -s, --size                 Specify the number of SRAs to use with Pilon,\n"
printf "                              this only applies when no SRAs are input\n"
printf "                              This will affect runtime, default is 5\n\n"
printf "   --outdir                   Specify output directory for files\n\n"
printf "   -k, --keep                 Keep the intermediate output files,\n"
printf "                              default is to remove files upon completion\n\n"

exit 1
}



while :; do
    casevar=`echo $1 | tr '[:upper:]' '[:lower:]'`
	case $casevar in
        -h | --h | --help | -help) display_help
                     exit 0
        ;;
        --genome | -genome) genome=$2
                 shift
                 param=$2
                 while [[ ${param:0:1} != '-' && ${param:0:1} != '' ]]; do
                     genome="$genome $2"
                     shift
                     param=$2
                 done
		;;
        --accession | -accession | --acc | -acc)   acc=true
                                                   echo $2 > "${DIR}/input_dna.out"
                                                   shift
                                                   param=$2
                                                   while [[ ${param:0:1} != '-' && ${param:0:1} != '' ]]; do
                                                        echo $2 >> "${DIR}/input_dna.out"
                                                        #input_dna="$input_dna $2"
                                                        shift
                                                        param=$2
                                                   done
		;;
		--organism | -organism) org=$2
                    shift
                    param=$2
                    while [[ ${param:0:1} != '-' && ${param:0:1} != '' ]]; do
                        org="$org $2"
                        shift
                        param=$2
                    done
		;;
        --outdir | -outdir) outdir=$2
        ;;
        --DNAfiles | -DNAfiles | --dnafiles | -dnafiles) dnafile=true
                                                        echo "DNA FILES"
                                                         echo $2 > "${DIR}/dna_files.out"
                                                         shift
                                                         param=$2
                                                         while [[ ${param:0:1} != '-' && ${param:0:1} != '' ]]; do
                                                            echo $2 >> "${DIR}/dna_files.out"
                                                            #input_dna="$input_dna $2"
                                                            shift
                                                            param=$2
                                                         done

        ;;
        --keep | -keep | --k | -k) keep=true
        ;;
        --size | -size | --s | -s) size=$2
        ;;
		*) break;
	esac
    shift
done


#check to make sure genome was specified
if [ -z "$genome" ]; then
    printf "\nERROR: No genome specified! Please use -h for help dialog.\n"
    exit 1
fi

#check to see if file exists
if [ ! -f $genome ]; then
    printf "\nGenome file could not be opened!\n"
    exit 1
fi

#if there is no outdir specified use current directory
#also check to see if directory specified exists -- make it if not
if [ -z "$outdir" ]; then
    outdir=$PWD
else
    if [ ! -d "$outdir" ]; then
        mkdir "$outdir"
    fi
fi


#if user wants to keep intermediate data, make a new directory to store outputs
if [ $keep == true ]; then
    mkdir "${outdir}/mod3" 2>/dev/null
    #need to move over inputdna and acc files if created
    if [ -e "${DIR}/input_dna.out" ]; then
        cp "${DIR}/input_dna.out" "${outdir}/mod3/"
    fi

    if [ -e "${DIR}/dna_files.out" ]; then
        cp "${DIR}/dna_files.out" "${outdir}/mod3/"
    fi

    DIR="${outdir}/mod3";

fi

echo "${DIR}/input_data.out"


#echo $genome
#echo $org

#make the SRA file either from searching (then) or from user input file (else)
#dont do any searching if files were input even when accession numbers weren't
if [ $acc == false ] && [ $dnafile == false ]; then

    #perform a search using Edirect tools
    esearch -db SRA -query "${org}[orgn] AND biomol dna[Properties] AND wgs[Assay Type]" | efetch -db SRA -format runinfo > "${DIR}/SRA_info.out"

    #check if the file is empty
    if ! [ -s "${DIR}/SRA_info.out" ]; then
            #this means either there is no SRA data or organism is spelled incorrectly
            printf "\nERROR: No SRA data found. Check input spelling!\n"
            exit 1
    fi

    #parse the SRA data file to pull out run ids
    awk -F "\"*,\"*" '{print $1}' "${DIR}/SRA_info.out" > "${DIR}/awkSRA1.out"

    #get rid of column headers and empty lines
    grep -v '^Run\|^$' "${DIR}/awkSRA1.out" > "${DIR}/awkSRA2.out"

    #we now want to only keep size specified amount of lines from SRAs
    SRA_len=`wc -l < "${DIR}/awkSRA2.out"`

    if [ $SRA_len > $size ]; then
        head -n $size "${DIR}/awkSRA2.out" > "${DIR}/awkSRA.out"
    else
        "${DIR}/awkSRA2.out" > "${DIR}/awkSRA.out"
    fi

    #now we have acc data so set to true
    acc=true
fi

if [ $acc == true ]; then
    cat "${DIR}/input_dna.out" > "${DIR}/awkSRA2.out"
    #we now want to only keep size specified amount of lines from SRAs
    SRA_len=`wc -l < "${DIR}/awkSRA2.out"`

    if [ $SRA_len > $size ]; then
        head -n $size "${DIR}/awkSRA2.out" > "${DIR}/awkSRA.out"
    else
        "${DIR}/awkSRA2.out" > "${DIR}/awkSRA.out"
    fi
fi



mkdir "${DIR}/blastDB" 2>/dev/null
#make a blast DB from the input genome
makeblastdb -in $genome -dbtype nucl -parse_seqids -out "${DIR}/blastDB/mod3_db.db" \
|| { printf '\nmakeblastdb failed! Please check input files and that magicblast is installed and in PATH variable!\n' ; exit 1; }

mkdir "${DIR}/SRA_aligned_mod1" 2>/dev/null

#loop through the SRAs line by line and align to genome
#we can let 4 run in parallel at once
if [ $acc == true ]; then
    N=4
    while read sra; do
        echo $sra
        ((i=i%N))
        ((i++==0)) && wait
        magicblast -db "${DIR}/blastDB/mod3_db.db" -sra $sra -splice F -no_unaligned -num_threads 4 -out "${DIR}/SRA_aligned_mod1/${sra}_align" \
        || { printf "\nSRA ${sra} failed! Please check input files and that magicblast is installed and in PATH variable!\nNote other SRAs may still be successfully running if more than one were specified.\n"; exit 1; } &
    done < "${DIR}/awkSRA.out"
    wait
fi

if [ $dnafile == true ]; then
    N=4
    while read dna; do
        echo $dna
        ((i=i%N))
        ((i++==0)) && wait
        magicblast -db "${DIR}/blastDB/mod3_db.db" -query $dna -splice F -no_unaligned -num_threads 4 -out "${DIR}/SRA_aligned_mod1/${dna}_align" -infmt fastq \
        || { printf "\nSRA $dna failed! Please check input files and that magicblast is installed and in PATH variable!\nNote other SRAs may still be successfully running if more than one were specified.\n" ; exit 1; } &
    done < "${DIR}/dna_files.out"
    wait
fi

#sort the BAM files and index them
mkdir "${DIR}/sortedBAM" 2>/dev/null

for f in "${DIR}/SRA_aligned_mod1/*" ; do
    outf=`echo $f | awk -F/ '{print $2}'`
    samtools sort $f -O BAM -o "${DIR}/sortedBAM/${outf}.sorted.bam" \
    || { printf '\nSAMtools failed! Please check that SAMtools is installed and in PATH variable!\n' ; exit 1; }

done

for f in "${DIR}/sortedBAM/*.bam" ; do
    samtools index $f \
    || { printf '\nSAMtools failed! Please check that SAMtools is installed and in PATH variable!\n' ; exit 1; }
done


#now we want to run Pilon with all of the SRA data

#need to make a string to pass to pilon of all of the SRA datafiles
frag_var=""
for f in "${DIR}/sortedBAM/*.bam" ; do
    frag_var="${frag_var} $f"
done


exit 1
#call Pilon
java -jar /software/pilon-1.22.jar --genome $genome --frags $frag_var --outdir $outdir \
|| { printf '\nPilon failed! Please check input and that Pilon is installed and in PATH variable!\n' ; exit 1; }

#exit
printf "\nThank you for using MASSIF BLAST\n"
exit 1



