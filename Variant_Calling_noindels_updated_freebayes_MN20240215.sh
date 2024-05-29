#!/bin/bash

echo -n "Complete path to run folder: "
read runfolder
echo -n "Complete path to raw data folder: "
read datafolder
echo -n "Reference file name (w/out .fasta extension): "
read reference
echo -n "Enter sample names in the format SAMPLENAME_SXX (separate samples with a space): "
read -a sample

echo "Select a mode for SNP calling (type 1-6): "
select SNP_Mode in 0.5_SNPs 0.2_SNPs 0.1_SNPs 0.05_SNPs 0.02_SNPs 0.01_SNPs Custom
do
	case $SNP_Mode in
		0.5_SNPs)
			mode="0.5_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.5" ;\
			;;
		0.2_SNPs)
			mode="0.2_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.2" ;\
			;;
		0.1_SNPs)
			mode="0.1_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.1" ;\
			;;
		0.05_SNPs)
			mode="0.05_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.05" ;\
			;;
		0.02_SNPs)
			mode="0.02_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.02" ;\
			;;
		0.01_SNPs)
			mode="0.01_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.01" ;\
			;;
		Custom)
			mode="Custom1" ;\
			echo "Specify the following parameters for variant calling: " ;\
			echo -n "Minimum mapping quality: " ;\
			read minmapq ;\
			echo -n "Minimum count of reads/observations to support an alternative allele: " ;\
			read mincount ;\
			echo -n "Minimum coverage for SNP calling: " ;\
			read mincov ;\
			echo -n "Minimum fraction of reads/observations to support an alternative allele (from 0 to 1): " ;\
			read minfrac ;\
			;;
		*) echo "Invalid Selection."
		;;
	esac
		for i in "${sample[@]}"
		do
		cd "$runfolder"/ ;\
		mkdir -p "$i"/{Raw_Data,QC,Trimming,Alignment,Variant_Calling_"$mode",Stats,"$i"_Export} ;\
		cp "$datafolder"/"$i"_L001_R1_001.fastq.gz "$datafolder"/"$i"_L001_R2_001.fastq.gz "$runfolder"/"$i"/ ;\
		cd ./"$i"/ ;\
		java -jar /home/fuberlin/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 "$i"_L001_R1_001.fastq.gz "$i"_L001_R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP://home/fuberlin/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;\
		mv "$runfolder"/"$i"/"$i"_forward_unpaired.fq.gz "$runfolder"/"$i"/"$i"_reverse_unpaired.fq.gz "$runfolder"/"$i"/Trimming/ ;\
		/home/fuberlin/Programs/fastqc_v0.11.9/FastQC/fastqc "$i"_forward_paired.fq.gz "$i"_reverse_paired.fq.gz --outdir="$runfolder"/"$i"/QC/ ;\
		rm "$runfolder"/"$i"/QC/"$i"_forward_paired_fastqc.zip "$runfolder"/"$i"/QC/"$i"_reverse_paired_fastqc.zip ;\
		mv "$i"_L001_R1_001.fastq.gz "$i"_L001_R2_001.fastq.gz "$runfolder"/"$i"/Raw_Data/ ;\
		mv "$i"_forward_paired.fq.gz "$i"_reverse_paired.fq.gz "$runfolder"/"$i"/Alignment/ ;\
		cd "$runfolder"/ ;\
		cp "$reference".fasta "$runfolder"/"$i"/Alignment/ ;\
		cd "$runfolder"/"$i"/Alignment/ ;\
		bwa index "$reference".fasta ;\
		time bwa mem "$reference".fasta "$i"_forward_paired.fq.gz "$i"_reverse_paired.fq.gz >bwa_paired_"$i".sam ;\
		samtools view -S -h -b bwa_paired_"$i".sam > bwa_paired_"$i".bam ;\
		samtools sort bwa_paired_"$i".bam -o bwa_paired_"$i"_sorted.bam ;\
		samtools index bwa_paired_"$i"_sorted.bam ;\
		cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/Variant_Calling_"$mode"/ ;\
		cd "$runfolder"/"$i"/Variant_Calling_"$mode"/ ;\
		cp "$runfolder"/"$i"/Alignment/"$reference".fasta "$runfolder"/"$i"/Variant_Calling_"$mode"/ ;\
		samtools faidx "$reference".fasta ;\
		conda init ;\
		conda activate freebayes ;\
		freebayes -f "$reference".fasta -m "$minmapq" -C "$mincount" --min-coverage "$mincov" -F "$minfrac" --pvar 0.001 -i bwa_paired_"$i"_sorted.bam >bwa_paired_"$i"_sorted_freebayes_"$mode".vcf ;\
		conda deactivate ;\
		cat bwa_paired_"$i"_sorted_freebayes_"$mode".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode".vcf.xls ;\
		bgzip -c bwa_paired_"$i"_sorted_freebayes_"$mode".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode".vcf.gz ;\
		tabix -f bwa_paired_"$i"_sorted_freebayes_"$mode".vcf.gz ;\
		bcftools consensus -f "$reference".fasta bwa_paired_"$i"_sorted_freebayes_"$mode".vcf.gz >consensus_"$i"_freebayes_"$mode".fa ;\
		samtools flagstat bwa_paired_"$i"_sorted.bam > "$i"_mapped.txt ;\
		echo "$(<"$i"_mapped.txt)" ;\
		mv "$i"_mapped.txt "$runfolder"/"$i"/Stats/ ;\
 		cd "$runfolder"/"$i"/Alignment/ ;\
		gunzip -k "$i"_forward_paired.fq.gz ;\
		wc -l "$i"_forward_paired.fq ;\
		cd "$runfolder"/"$i"/Variant_Calling_"$mode"/ ;\
		samtools coverage -m bwa_paired_"$i"_sorted.bam -o coverage_"$i".txt ;\
		samtools depth bwa_paired_"$i"_sorted.bam | awk '{sum+=$3} END {print "Average = ", sum/NR}' > avgdepth_"$i".txt ;\
		samtools depth -a bwa_paired_"$i"_sorted.bam > depth_"$i".csv ;\
		cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/"$i"_Export/ ;\
		cp "$runfolder"/"$i"/Alignment/bwa_paired_"$i"_sorted.bam.bai "$runfolder"/"$i"/"$i"_Export/ ;\
		cp bwa_paired_"$i"_sorted_freebayes_"$mode".vcf "$runfolder"/"$i"/"$i"_Export/ ;\
		cp bwa_paired_"$i"_sorted_freebayes_"$mode".vcf.xls "$runfolder"/"$i"/"$i"_Export/ ;\
		cp consensus_"$i"_freebayes_"$mode".fa "$runfolder"/"$i"/"$i"_Export/ ;\
		mv depth_"$i".csv "$runfolder"/"$i"/Stats/ ;\
		mv avgdepth_"$i".txt "$runfolder"/"$i"/Stats/ ;\
		mv coverage_"$i".txt "$runfolder"/"$i"/Stats/ ;\
		cd "$runfolder"/"$i"/Stats/ ;\
		echo "$(<avgdepth_"$i".txt)" ;\
		cd "$runfolder"/"$i"/ ;\
		cp -r "$runfolder"/"$i"/Stats/ "$runfolder"/"$i"/"$i"_Export/ ;\
		done
	echo "Re-do Variant Calling with different parameters? (type 1-7): " ;\
	select SNP_Mode2 in No 0.5_SNPs 0.2_SNPs 0.1_SNPs 0.05_SNPs 0.02_SNPs 0.01_SNPs Custom ;\
	do
	case $SNP_Mode2 in
		No)
			for i in "${sample[@]}"
			do
				rm -r "$runfolder"/"$i"/Raw_Data/ ;\
				rm -r "$runfolder"/"$i"/QC/ ;\
				rm -r "$runfolder"/"$i"/Trimming/ ;\
				rm -r "$runfolder"/"$i"/Alignment/ ;\
				rm -r "$runfolder"/"$i"/Stats/ ;\
			done
			exit 0
			;;
		0.5_SNPs)
			mode2="0.5_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.5" ;\
			;;
		0.2_SNPs)
			mode2="0.2_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.2" ;\
			;;
		0.1_SNPs)
			mode2="0.1_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.1" ;\
			;;
		0.05_SNPs)
			mode2="0.05_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.05" ;\
			;;
		0.02_SNPs)
			mode2="0.02_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.02" ;\
			;;
		0.01_SNPs)
			mode2="0.01_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.01" ;\
			;;
		Custom)
			mode2="Custom2" ;\
			echo "Specify the following parameters for variant calling: " ;\
			echo -n "Minimum mapping quality: " ;\
			read minmapq ;\
			echo -n "Minimum count of reads/observations to support an alternative allele: " ;\
			read mincount ;\
			echo -n "Minimum coverage for SNP calling: " ;\
			read mincov ;\
			echo -n "Minimum fraction of reads/observations to support an alternative allele (from 0 to 1): " ;\
			read minfrac ;\
			;;
		*) echo "Invalid Selection."
		;;
	esac
		for i in "${sample[@]}"
		do
			cd "$runfolder"/"$i"/ ;\
			mkdir -p Variant_Calling_"$mode2" ;\
			cd "$runfolder"/"$i"/Alignment/ ;\
			cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/Variant_Calling_"$mode2"/ ;\
			cd "$runfolder"/"$i"/Variant_Calling_"$mode2"/ ;\
			cp "$runfolder"/"$i"/Alignment/"$reference".fasta "$runfolder"/"$i"/Variant_Calling_"$mode2"/ ;\
			samtools faidx "$reference".fasta ;\
			conda init ;\
			conda activate freebayes ;\
			freebayes -f "$reference".fasta -m "$minmapq" -C "$mincount" --min-coverage "$mincov" -F "$minfrac" --pvar 0.001 -i bwa_paired_"$i"_sorted.bam >bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf ;\
			conda deactivate ;\
			cat bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf.xls ;\
			bgzip -c bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf.gz ;\
			tabix -f bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf.gz ;\
			bcftools consensus -f "$reference".fasta bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf.gz >consensus_"$i"_freebayes_"$mode2".fa ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf "$runfolder"/"$i"/"$i"_Export/ ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode2".vcf.xls "$runfolder"/"$i"/"$i"_Export/ ;\
			cp consensus_"$i"_freebayes_"$mode2".fa "$runfolder"/"$i"/"$i"_Export/ ;\
		done
	done
	echo "Re-do Variant Calling with different parameters? (type 1-7): " ;\
	select SNP_Mode3 in No 0.5_SNPs 0.2_SNPs 0.1_SNPs 0.05_SNPs 0.02_SNPs 0.01_SNPs Custom ;\
	do
	case $SNP_Mode3 in
		No)
			for i in "${sample[@]}"
			do
				rm -r "$runfolder"/"$i"/Raw_Data/ ;\
				rm -r "$runfolder"/"$i"/QC/ ;\
				rm -r "$runfolder"/"$i"/Trimming/ ;\
				rm -r "$runfolder"/"$i"/Alignment/ ;\
				rm -r "$runfolder"/"$i"/Stats/ ;\
			done
			exit 0
			;;
		0.5_SNPs)
			mode3="0.5_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.5" ;\
			;;
		0.2_SNPs)
			mode3="0.2_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.2" ;\
			;;
		0.1_SNPs)
			mode3="0.1_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.1" ;\
			;;
		0.05_SNPs)
			mode3="0.05_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.05" ;\
			;;
		0.02_SNPs)
			mode3="0.02_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.02" ;\
			;;
		0.01_SNPs)
			mode3="0.01_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.01" ;\
			;;
		Custom)
			mode3="Custom3" ;\
			echo "Specify the following parameters for variant calling: " ;\
			echo -n "Minimum mapping quality: " ;\
			read minmapq ;\
			echo -n "Minimum count of reads/observations to support an alternative allele: " ;\
			read mincount ;\
			echo -n "Minimum coverage for SNP calling: " ;\
			read mincov ;\
			echo -n "Minimum fraction of reads/observations to support an alternative allele (from 0 to 1): " ;\
			read minfrac ;\
			;;
		*) echo "Invalid Selection."
		;;
	esac
		for i in "${sample[@]}"
		do
			cd "$runfolder"/"$i"/ ;\
			mkdir -p Variant_Calling_"$mode3" ;\
			cd "$runfolder"/"$i"/Alignment/ ;\
			cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/Variant_Calling_"$mode3"/ ;\
			cd "$runfolder"/"$i"/Variant_Calling_"$mode3"/ ;\
			cp "$runfolder"/"$i"/Alignment/"$reference".fasta "$runfolder"/"$i"/Variant_Calling_"$mode3"/ ;\
			samtools faidx "$reference".fasta ;\
			conda init ;\
			conda activate freebayes ;\
			freebayes -f "$reference".fasta -m "$minmapq" -C "$mincount" --min-coverage "$mincov" -F "$minfrac" --pvar 0.001 -i bwa_paired_"$i"_sorted.bam >bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf ;\
			conda deactivate ;\
			cat bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf.xls ;\
			bgzip -c bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf.gz ;\
			tabix -f bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf.gz ;\
			bcftools consensus -f "$reference".fasta bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf.gz >consensus_"$i"_freebayes_"$mode3".fa ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf "$runfolder"/"$i"/"$i"_Export/ ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode3".vcf.xls "$runfolder"/"$i"/"$i"_Export/ ;\
			cp consensus_"$i"_freebayes_"$mode3".fa "$runfolder"/"$i"/"$i"_Export/ ;\
		done
	done
	echo "Re-do Variant Calling with different parameters? (type 1-7): " ;\
	select SNP_Mode4 in No 0.5_SNPs 0.2_SNPs 0.1_SNPs 0.05_SNPs 0.02_SNPs 0.01_SNPs Custom ;\
	do
	case $SNP_Mode4 in
		No)
			for i in "${sample[@]}"
			do
				rm -r "$runfolder"/"$i"/Raw_Data/ ;\
				rm -r "$runfolder"/"$i"/QC/ ;\
				rm -r "$runfolder"/"$i"/Trimming/ ;\
				rm -r "$runfolder"/"$i"/Alignment/ ;\
				rm -r "$runfolder"/"$i"/Stats/ ;\
			done
			exit 0
			;;
		0.5_SNPs)
			mode4="0.5_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.5" ;\
			;;
		0.2_SNPs)
			mode4="0.2_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.2" ;\
			;;
		0.1_SNPs)
			mode4="0.1_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.1" ;\
			;;
		0.05_SNPs)
			mode4="0.05_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.05" ;\
			;;
		0.02_SNPs)
			mode4="0.02_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.02" ;\
			;;
		0.01_SNPs)
			mode4="0.01_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.01" ;\
			;;
		Custom)
			mode4="Custom4" ;\
			echo "Specify the following parameters for variant calling: " ;\
			echo -n "Minimum mapping quality: " ;\
			read minmapq ;\
			echo -n "Minimum count of reads/observations to support an alternative allele: " ;\
			read mincount ;\
			echo -n "Minimum coverage for SNP calling: " ;\
			read mincov ;\
			echo -n "Minimum fraction of reads/observations to support an alternative allele (from 0 to 1): " ;\
			read minfrac ;\
			;;
		*) echo "Invalid Selection."
		;;
	esac
		for i in "${sample[@]}"
		do
			cd "$runfolder"/"$i"/ ;\
			mkdir -p Variant_Calling_"$mode4" ;\
			cd "$runfolder"/"$i"/Alignment/ ;\
			cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/Variant_Calling_"$mode4"/ ;\
			cd "$runfolder"/"$i"/Variant_Calling_"$mode4"/ ;\
			cp "$runfolder"/"$i"/Alignment/"$reference".fasta "$runfolder"/"$i"/Variant_Calling_"$mode4"/ ;\
			samtools faidx "$reference".fasta ;\
			conda init ;\
			conda activate freebayes ;\
			freebayes -f "$reference".fasta -m "$minmapq" -C "$mincount" --min-coverage "$mincov" -F "$minfrac" --pvar 0.001 -i bwa_paired_"$i"_sorted.bam >bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf ;\
			conda deactivate ;\
			cat bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf.xls ;\
			bgzip -c bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf.gz ;\
			tabix -f bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf.gz ;\
			bcftools consensus -f "$reference".fasta bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf.gz >consensus_"$i"_freebayes_"$mode4".fa ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf "$runfolder"/"$i"/"$i"_Export/ ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode4".vcf.xls "$runfolder"/"$i"/"$i"_Export/ ;\
			cp consensus_"$i"_freebayes_"$mode4".fa "$runfolder"/"$i"/"$i"_Export/ ;\
		done
	done
	echo "Re-do Variant Calling with different parameters? (type 1-7): " ;\
	select SNP_Mode5 in No 0.5_SNPs 0.1_SNPs 0.2_SNPs 0.05_SNPs 0.02_SNPs 0.01_SNPs Custom ;\
	do
	case $SNP_Mode5 in
		No)
			for i in "${sample[@]}"
			do
				rm -r "$runfolder"/"$i"/Raw_Data/ ;\
				rm -r "$runfolder"/"$i"/QC/ ;\
				rm -r "$runfolder"/"$i"/Trimming/ ;\
				rm -r "$runfolder"/"$i"/Alignment/ ;\
				rm -r "$runfolder"/"$i"/Stats/ ;\
			done
			exit 0
			;;
		0.5_SNPs)
			mode5="0.5_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.5" ;\
			;;
		0.2_SNPs)
			mode5="0.2_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.2" ;\
			;;
		0.1_SNPs)
			mode5="0.1_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.1" ;\
			;;
		0.05_SNPs)
			mode5="0.05_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.05" ;\
			;;
		0.02_SNPs)
			mode5="0.02_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.02" ;\
			;;
		0.01_SNPs)
			mode5="0.01_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.01" ;\
			;;
		Custom)
			mode5="Custom5" ;\
			echo "Specify the following parameters for variant calling: " ;\
			echo -n "Minimum mapping quality: " ;\
			read minmapq ;\
			echo -n "Minimum count of reads/observations to support an alternative allele: " ;\
			read mincount ;\
			echo -n "Minimum coverage for SNP calling: " ;\
			read mincov ;\
			echo -n "Minimum fraction of reads/observations to support an alternative allele (from 0 to 1): " ;\
			read minfrac ;\
			;;
		*) echo "Invalid Selection."
		;;
	esac
		for i in "${sample[@]}"
		do
			cd "$runfolder"/"$i"/ ;\
			mkdir -p Variant_Calling_"$mode5" ;\
			cd "$runfolder"/"$i"/Alignment/ ;\
			cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/Variant_Calling_"$mode5"/ ;\
			cd "$runfolder"/"$i"/Variant_Calling_"$mode5"/ ;\
			cp "$runfolder"/"$i"/Alignment/"$reference".fasta "$runfolder"/"$i"/Variant_Calling_"$mode5"/ ;\
			samtools faidx "$reference".fasta ;\
			conda init ;\
			conda activate freebayes ;\
			freebayes -f "$reference".fasta -m "$minmapq" -C "$mincount" --min-coverage "$mincov" -F "$minfrac" --pvar 0.001 -i bwa_paired_"$i"_sorted.bam >bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf ;\
			conda deactivate ;\
			cat bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf.xls ;\
			bgzip -c bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf.gz ;\
			tabix -f bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf.gz ;\
			bcftools consensus -f "$reference".fasta bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf.gz >consensus_"$i"_freebayes_"$mode5".fa ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf "$runfolder"/"$i"/"$i"_Export/ ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode5".vcf.xls "$runfolder"/"$i"/"$i"_Export/ ;\
			cp consensus_"$i"_freebayes_"$mode5".fa "$runfolder"/"$i"/"$i"_Export/ ;\
		done
	done
	echo "Re-do Variant Calling with different parameters? (type 1-7): " ;\
	select SNP_Mode6 in No 0.5_SNPs 0.2_SNPs 0.1_SNPs 0.05_SNPs 0.02_SNPs 0.01_SNPs Custom ;\
	do
	case $SNP_Mode6 in
		No)
			for i in "${sample[@]}"
			do
				rm -r "$runfolder"/"$i"/Raw_Data/ ;\
				rm -r "$runfolder"/"$i"/QC/ ;\
				rm -r "$runfolder"/"$i"/Trimming/ ;\
				rm -r "$runfolder"/"$i"/Alignment/ ;\
				rm -r "$runfolder"/"$i"/Stats/ ;\
			done
			exit 0
			;;
		0.5_SNPs)
			mode6="0.5_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.5" ;\
			;;
		0.2_SNPs)
			mode6="0.2_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.2" ;\
			;;
		0.1_SNPs)
			mode6="0.1_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.1" ;\
			;;
		0.05_SNPs)
			mode6="0.05_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.05" ;\
			;;
		0.02_SNPs)
			mode6="0.02_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.02" ;\
			;;
		0.01_SNPs)
			mode6="0.01_SNPs" ;\
			minmapq="5" ;\
			mincount="3" ;\
			mincov="10" ;\
			minfrac="0.01" ;\
			;;
		Custom)
			mode6="Custom5" ;\
			echo "Specify the following parameters for variant calling: " ;\
			echo -n "Minimum mapping quality: " ;\
			read minmapq ;\
			echo -n "Minimum count of reads/observations to support an alternative allele: " ;\
			read mincount ;\
			echo -n "Minimum coverage for SNP calling: " ;\
			read mincov ;\
			echo -n "Minimum fraction of reads/observations to support an alternative allele (from 0 to 1): " ;\
			read minfrac ;\
			;;
		*) echo "Invalid Selection."
		;;
	esac
		for i in "${sample[@]}"
		do
			cd "$runfolder"/"$i"/ ;\
			mkdir -p Variant_Calling_"$mode6" ;\
			cd "$runfolder"/"$i"/Alignment/ ;\
			cp bwa_paired_"$i"_sorted.bam "$runfolder"/"$i"/Variant_Calling_"$mode6"/ ;\
			cd "$runfolder"/"$i"/Variant_Calling_"$mode6"/ ;\
			cp "$runfolder"/"$i"/Alignment/"$reference".fasta "$runfolder"/"$i"/Variant_Calling_"$mode6"/ ;\
			samtools faidx "$reference".fasta ;\
			conda init ;\
			conda activate freebayes ;\
			freebayes -f "$reference".fasta -m "$minmapq" -C "$mincount" --min-coverage "$mincov" -F "$minfrac" --pvar 0.001 -i bwa_paired_"$i"_sorted.bam >bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf ;\
			conda deactivate ;\
			cat bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf.xls ;\
			bgzip -c bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf >bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf.gz ;\
			tabix -f bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf.gz ;\
			bcftools consensus -f "$reference".fasta bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf.gz >consensus_"$i"_freebayes_"$mode6".fa ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf "$runfolder"/"$i"/"$i"_Export/ ;\
			cp bwa_paired_"$i"_sorted_freebayes_"$mode6".vcf.xls "$runfolder"/"$i"/"$i"_Export/ ;\
			cp consensus_"$i"_freebayes_"$mode6".fa "$runfolder"/"$i"/"$i"_Export/ ;\
		done
		exit 0
	done 
	exit 0
done
exit 0
