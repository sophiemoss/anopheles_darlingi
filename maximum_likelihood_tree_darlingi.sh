# create reference sequence for mito genome
# samtools faidx AnoDarl_H01.genomic.fasta NC_064612.1 > AnoDarl_H01_NC_064612.1.fasta
# samtools faidx AnoDarl_H01_NC_064612.1.fasta
#  use /mnt/storage11/sophie/gitrepos/anopheles_darlingi/create_mito_fasta_sequences.sh
# this will generate a gvcf folder with the sequences for each sample, concatenate all of these together.
# cat *_consensus.fasta > darlingi_mito_consensus_sequences.fasta
# download and view the alignment in aliview, and align

# Both types of tree
raxml-ng --all --msa darlingi_allsamples_consensus_mito.aligned.fasta --model GTR --prefix Anopheles_darlingi --seed 729264 --bs-metric tbe --tree rand{1} --bs-trees 1000

# Use best_tree for iTOL
