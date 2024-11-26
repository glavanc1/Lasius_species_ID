# Competitive Mapping

## Concatenate target and contaminant genomes
` cat *.fasta > cat_gen.fasta `

## Map individuals to concatenated genome
` bwa index cat_gen.fasta `  
` bwa mem -t 20 cat_gen.fasta $READ.fq.gz| samtools sort --output-fmt=SAM --threads 20 | samtools view --no-header -o ./sam/catgen_$READ.sam - `  

## Identify reads that did map to *Lasius*
` Rscript filtercontaminates.r `

## Filter out reads that did not map to *Lasius*
` filterbyname.sh in=${SAMPS}/${READ}.fq.gz names=${IDDIR}/${READ}_ID.txt out=${OUT}/${READ}.filtered.fq.gz include=t `
