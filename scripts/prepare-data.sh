
. scripts/setup.sh

# Create output folder
mkdir -p refs

# VectorBase - Aedes Gene Feature File GTF
echo "*** [$(date)] [prepare-data.sh] GTF Setup"
if [ ! -f $GTF ]; then
  echo "*** [$(date)] Download GTF"
  curl -L https://www.vectorbase.org/download/aedes-aegypti-lvpagwgbasefeaturesaaegl51gtfgz -o refs/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1.gtf.gz
  gunzip ${GTF}.gz
fi

# VectorBase - Aedes Reference Genome
echo "*** [$(date)] [prepare-data.sh] Genome Setup"
if [ ! -f $REF ]; then
  echo "*** [$(date)] Download Genome"
  curl -L https://www.vectorbase.org/download/aedes-aegypti-lvpagwgchromosomesaaegl5fagz -o refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz
  gunzip ${REF}.gz

fi

echo "*** [$(date)] [prepare-data.sh] Index Genome - $ALIGNER"
if [ $ALIGNER = "star" ]; then
  STAR --runThreadN 8 --runMode genomeGenerate  \
  --genomeDir /work/refs \
  --genomeFastaFiles $REF \
  --sjdbGTFfile $GTF \
  --sjdbOverhang 48 \
  --limitGenomeGenerateRAM 28000000000
else
  hisat2-build $REF $REF
fi


###### TODO -> Download Sample Data

# Unzip Sample data
echo "*** [$(date)] [prepare-data.sh] Unzip sample data"
SAMPLES=($(ls reads/*.fastq.gz 2> /dev/null))
for SAMPLE in ${SAMPLES[@]};
do
    gunzip $SAMPLE
done

echo "*** [$(date)] [prepare-data.sh] Done"
