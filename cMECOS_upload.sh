rawfq1=./example/rawdata/sample/V350110680_L01_1_1.fq.gz
rawfq2=./example/rawdata/sample/V350110680_L01_1_2.fq.gz
sample=sample

###step6.contig binning
cd ./example/03.binning/$sample
ln -s ./example/02.assembly/$sample/athena/results/olc/athena.asm.fa  $sample.fasta
bowtie2-build  -f $sample.fasta $sample.fasta
bowtie2 -x $sample.fasta  -1 ./example/02.assembly/$sample/$sample.unmap.clean.1.fq.paired.fq  -2 ./example/02.assembly/$sample/$sample.unmap.clean.2.fq.paired.fq  -p 20 -S $sample.sam 2> bowtie2.log
samtools view -F 4  -Sb $sample.sam  > $sample.bam
samtools  sort $sample.bam  -o $sample.sort.bam
samtools index  $sample.sort.bam
jgi_summarize_bam_contig_depths  --outputDepth $sample.depth.txt $sample.sort.bam
metabat2 -i $sample.fasta -a $sample.depth.txt -o $sample --sensitive -t 20 -v > $sample.log.txt

###step7.prokka annotation
ls ./example/03.binning/$sample/$sample.*.fa >./example/04.annotation/$sample.binning.list
for $i in $sample.binning.list
do
prokka --force --cpus 10 --outdir ./example/04.annotation/prokka/$sample/prokka/$i --prefix $i  --locustag $i  --metagenome --kingdom Bacteria  ./example/03.binning/$sample/$i
done

###step8.emapper annotation
for $i in $sample.binning.list
do
emapper.py -i $i --itype metagenome -m diamond --evalue 1e-05 -o $i --output_dir ./example/04.annotation/emapper/$sample/output
done

###step9.gtdb annotation
ln -s ./example/03.binning/$sample/$sample.*.fa  ./example/04.annotation/gtdb/contig/
gtdbtk  classify_wf --genome_dir ./example/04.annotation/gtdb/contig/   --out_dir ./example/04.annotation/gtdb/new_gtdb --extension fa --cpus 30

###step10.vfdb annotation
for $i in $sample.binning.list
do
abricate $i --db vfdb --minid=75 > ./example/04.annotation/vfdb/$sample/$i.vfdb.tab
done

###step11.resfinder annotation
for $i in $sample.binning.list
do
abricate $i --db resfinder_new --minid=75 > ./example/04.annotation/resfinder/$sample/$i.resfinder_new.tab
done
