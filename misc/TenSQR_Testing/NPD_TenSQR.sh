#!/bin/sh
#SBATCH --mem=40g

python=/export/home/jhe/download/anaconda/anaconda3/envs/tensqr_env/bin/python
ExtractMatrix=/export/home/jhe/download/TenSQR/TenSQR-master/ExtractMatrix
java=/export/home/jhe/download/java_jdk_8u201/jdk1.8.0_201/bin/java
clustalo=/export/home/jhe/.local/bin/clustalo

inputdir=$1
inputdir+="/input/"
interdir=$1
interdir+="/intermediate/"
gsdir=$1
gsdir+="/gold_standard/"
outdir=$1
outdir+="/output/"
prefix=$2
seqerr=$3
popsize=$4
pools=$5

prefixhead="$interdir$prefix"
refseq="$inputdir"
refseq+="HIV_HXB2"

ownpdir="/gpfs/home/lmak/PHX_NonPerfect_Data4/programs/"

for (( p=0; p<$pools; p++));
do
	confighead=$prefixhead
	confighead+="_$p"
	samhead=$prefixhead
	samhead+="_p$p"
	/gpfs/common/programs/samtools-0.1.19 view -ho "$samhead".sam "$samhead".bam
	sed -e "s|rep1|$samhead|g;s|rep2|$seqerr|g;s|rep3|$popsize|g" "$inputdir"TenSQR.properties > "$confighead"_TenSQR.properties
	$ExtractMatrix "$confighead"_TenSQR.properties
	$python "$ownpdir"TenSQR.py "$confighead"_TenSQR.properties
	$java -jar "$ownpdir"T2FA.jar "$refseq" "$samhead"
	$clustalo -i "$samhead"_ViralSeq.fasta -o "$samhead".fa
	$java -jar "$ownpdir"TSR2SD.jar "$samhead" 
done

# NOTE: In-pool inter_freq_vars will be located in intermediate/ because of how JH's program is set up. 
# TSRMerge.jar will put the global final inter_freq_vars and intra_freq in output/. O2R3.jar thus still output to output/.
$java -jar "$ownpdir"TSRMerge.jar $interdir $prefix $outdir $pools
$java -jar "$ownpdir"O2R3.jar $gsdir $outdir $prefix $pools 0.05
$java -jar "$ownpdir"O2R3.jar $gsdir $outdir $prefix $pools 0.10