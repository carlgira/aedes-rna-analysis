
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

  echo "*** [$(date)] [prepare-data.sh] Index Genome - $ALIGNER"
  if [ $ALIGNER = "star" ]; then
    STAR --runThreadN $NCPU --runMode genomeGenerate  \
    --genomeDir /work/refs \
    --genomeFastaFiles $REF \
    --sjdbGTFfile $GTF \
    --sjdbOverhang 48 \
    --limitGenomeGenerateRAM 26000000000
  else
    hisat2-build $REF $REF
  fi

fi

###### TODO -> Download Sample Data

# Unzip Sample data
echo "*** [$(date)] [prepare-data.sh] Unzip sample data"
SAMPLES=($(ls /work/reads/*.fq.gz 2> /dev/null))
echo "$SAMPLES"
for SAMPLE in ${SAMPLES[@]};
do
    gunzip $SAMPLE
done

echo "*** [$(date)] [prepare-data.sh] Done"
