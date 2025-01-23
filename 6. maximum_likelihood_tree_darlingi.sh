# create reference sequence for mito genome
# samtools faidx AnoDarl_H01.genomic.fasta NC_064612.1 > AnoDarl_H01_NC_064612.1.fasta
# samtools faidx AnoDarl_H01_NC_064612.1.fasta
#  use /mnt/storage11/sophie/gitrepos/anopheles_darlingi/create_mito_fasta_sequences.sh
# this will generate a gvcf folder with the sequences for each sample, and then concatenate all of these together.
# download and view the alignment in aliview, and align. Or align in terminal.

# large sequences so used mauve


# Both types of tree

# Mito only
raxml-ng --all --msa darlingi_phenotyped_samples_consensus_mitochondria.aligned.fasta  --model GTR --prefix Anopheles_darlingi_phenotyped --seed 729264 --bs-metric tbe --tree rand{1} --bs-trees 1000

# Whole genome
raxml-ng --all --msa allsamples_consensus_wholegenome.fasta --model GTR --prefix Anopheles_darlingi_wholegenome --seed 827435 --bs-metric tbe --tree rand{1} --bs-trees 1000


# Use best_tree for iTOL
