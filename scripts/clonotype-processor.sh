#!/bin/bash


sra=""
name=""
dir="default/path/to/copy/clonotypes/to"
cores=8
mem=30
startPath=$(pwd)

while getopts "s:n:hd:" args; do
	case "${args}" in
		s)
			sra=${OPTARG}
			;;
		n)
			name=${OPTARG}
			;;
		d)
			dir=${OPTARG}
			;;
		c)
			cores=${OPTARG}
			;;
		m)
			mem=${OPTARG}
			;;
		h)
			echo -e "This script initialises a pipeline in which SRA experiment with certain ID is downloaded, and processed with cellranger vdj and harvests the clonotypes from said analysis and puts them in a specified folder, if no directory is specified, the default is used. The default path is put in \"dir\" variable. It is necessery to pass both SRA ID and the name under which the clonotypes file will be saved so that one can make distintion between them.\n"
			echo -e "\t-s\n\t\tSRA ID. It is a mandatory argument.\n"
			echo -e "\t-n\n\t\tName of the clonotype file, can be the ailment of patient from whom the sample was taken. It is a mandatory argument.\n"
			echo -e "\t-d\n\t\tDirectory to which clonotypes should be moved. Optional argument.\n"
			echo -e "\t-c\n\t\tNumber of cores used for processing. Optional argument.\n"
			echo -e "\t-m\n\t\tAmount of memory, in gigabytes, used for processing. Optional argument.\n"
			echo -e "\t-h\n\t\tView this information.\n"
			exit 1
			;;
		:)
			echo -e "Error: -${OPTARG} requires an argument."
			exit 1                       
			;;
	esac
done

if [ $name = "" ] || [ $sra = "" ] 
then
	echo "Error: not enough arguments were provided"
	exit 1
fi

mkdir cellranger-working-dir
cd cellranger-working-dir
mkdir dataset-vdj
cd dataset-vdj
fastq-dump --split-files -v -L info $sra
for file in $(ls)
do
	gzip $file
	mv ${file}.gz ${sra}_S1_L001_R${file:12:1}_001.fastq.gz
done
cd ..
curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
cellranger vdj --id=$name --reference=refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 --fastqs=dataset-vdj --sample=$sra --localcores=8 --localmem=30
cd ${name}/outs
mv clonotypes.csv ${name}_clonotypes.csv
cp ${name}_clonotypes.csv ${dir}
cd $startPath
rm -r cellranger-working-dir
