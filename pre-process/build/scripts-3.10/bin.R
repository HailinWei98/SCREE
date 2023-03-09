library(Signac)

args <- commandArgs(T)
fragments <- args[1]
binsize <- args[2]
process_n <- args[3]
genome <- args[4]
file.path <- args[5]


fragments <- CreateFragmentObject(fragments)
genome <- read.table(genome, header = F)
genome_name <- genome[, 1]
genome <- as.vector(genome[, 2])
names(genome) <- genome_name

a <- GenomeBinMatrix(
    fragments = fragments,
    genome = genome,
    binsize = as.numeric(binsize),
    process_n = as.numeric(process_n),
)

saveRDS(a, file = file.path)