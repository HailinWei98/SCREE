import os
version = "0.0.4"
def main():
    #from optparse import OptionParser
    import argparse
    #from pkg_resources import resource_filename
    import sys
    #para = OptionParser(usage = 'SCREE [--help/-h][--datatype][--filetype][--reference][--config][--mkref][--input][--sample][--output][-m][-c][--csv]')
    parser = argparse.ArgumentParser(description = 'SCREE(Single-cell CRISPR scREen data analyses and pErturbation modeling) is a workflow for single-cell CRISPR screens data processing and analysis.')
    parser.add_argument("-v", "--version", action = "store_true", help = "Print version info.")
    subparsers = parser.add_subparsers(dest = "subcommand")
    scree_preprocess = subparsers.add_parser('preprocess', help = 'Reads alignment and quantification.')
    scree_analysis = subparsers.add_parser('analysis', help = 'Performing all downstream analyses with one command.')
    scree_preprocess.add_argument("--datatype", default = "RNA", choices = ["RNA", "ATAC"], type = str, help = "Input data type, can be one of 'RNA', 'ATAC', corresponding to scRNA-seq and scATAC-seq, respectively.")
    scree_preprocess.add_argument("--filetype", default = "multi", choices = ["multi", "single"], type = str, help = "Input file type, one of 'multi', 'single'. 'multi' indicates sgRNA sequence and mRNA/DNA sequence are not in the same files, 'single' indicates sgRNA sequence and mRNA sequence are in the same files. Notably, for scATAC-seq input, SCREE only support 'multi'.")
#    para.add_option("--fasta", default = "", type = str, help = "File path of reference in FASTA format.")
#    para.add_option("--gtf", "-g", default = "RNA", type = str, help = "Input data type, can be one of #RNA/ATAC, only support scATAC-seq or scRNA-seq.")
    scree_preprocess.add_argument("--reference", default = "", dest = "reference", type = str, help = "Reference path.")
    scree_preprocess.add_argument("--config", default = "", dest = "config", type = str, help = "Config file to generate reference index.")
    scree_preprocess.add_argument("--mkref", action = "store_true", dest = "mkref", help = "Generate reference index or not.")
    scree_preprocess.add_argument("--input", default = "", dest = "input", type = str, help = "File path of all the fastq files.")
    scree_preprocess.add_argument("--sample", default = "", dest = "sample", type = str, help = "Sample prefix.")
    scree_preprocess.add_argument("--output", default = ".", dest = "output", type = str, help = "Output folder to create.")
    scree_preprocess.add_argument("-m", dest = "memory", default = 50, type = int, help = "Maximum local memory to use. Unit: GiB")
    scree_preprocess.add_argument("-c", dest = "cores", default = 8, type = int, help = "Maximum local cores to use. ")
    scree_preprocess.add_argument("--csv", default = "", dest = "csv", type = str, help = "Config file to process 'multi' file type.")
    scree_preprocess.add_argument("--bin", action = "store_true", dest = "bin", help = "Generate a bin-based matrix instead of a peak-based matrix.")
    scree_preprocess.add_argument("--binsize", dest = "binsize", default = 5000, type = int, help = "Size of the genome bins to use. Only available with the parameter '--bin'. Default is 5000. ")
    scree_preprocess.add_argument("-n", "--process_n", dest = "process_n", default = 2000, type = int, help = "Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory. Only available with the parameter '--bin'. Default is 2000. ")
    scree_preprocess.add_argument("-l", "--chrLength", dest = "chrLength", type = str, help = "Path to the file which stores the chromosome length of the reference genome, always named as 'chrNameLength.txt'. Only available with the parameter '--bin'.")
    scree_preprocess.add_argument("--para", default = "", dest = "para", type = str, help = "Other parameters passed into cellranger or cellranger-atac, should be in string format like '--parameter1 value1 --parameter2 value2'. ")
    
    scree_analysis.add_argument("--type", dest = "type", type = str, default = "RNA", choices = ["RNA", "ATAC", "enhancer"], help = "Data type, one of 'RNA', 'ATAC', 'enhancer'. Default is 'RNA'.")
    scree_analysis.add_argument("--replicate", dest = "replicate", default = 1, help = "Directory to a txt file only contains the replicate information of each cell, in the same order of cells in the input matrix. If no replicate information, we will consider that all cells are from the same replicate, and this parameter will be set as 1. Default is 1.")
    scree_analysis.add_argument("--mtx", dest = "mtx", type = str, help = "Directory to an rds file of the SeuratObject, with cell in columns and features in rows, or the directory to the folder of 10X-like outputs.")
    scree_analysis.add_argument("--sg_lib", dest = "sg_lib", required = False, type = str, help = "Directory to a txt file of the sgRNA information. This parameter can be ommited when the 'sg_format' is '10X'. ")
    scree_analysis.add_argument("--sg_format", dest = "sg_format", type = str, default = "10X", choices = ["10X", "dataframe", "matrix", "table"], help = "Format of the input sg_lib, one of '10X', 'dataframe', 'matrix', 'table'. '10X' means the sgRNA information is an sgRNA count matrix stored in the 10X-like outputs; 'dataframe' means a table with 3 columns: cell(cell barcode), barcode(sgRNA name), freqUMI counts in each cell for each sgRNA; 'matrix' means a sgRNA count matrix with cells in columns and sgRNAs in rows; 'table' means a table with 3 columns: cell(cell barcode), barcode(sgRNA name), gene(sgRNA targeted gene). Default is '10X'. ")
    scree_analysis.add_argument("--freq_cut", dest = "freq_cut", default = 20, type = int, help = "Cutoff of sgRNA count numbers. For each cell, only sgRNA with counts more than freq_cut will be retained. Default is 20. ")
    scree_analysis.add_argument("--freq_percent", dest = "freq_percent", default = 0.8, type = float, help = "Cutoff of sgRNA frequency. For each cell, only sgRNA with frequency more than freq_percent will be retained. Default is 0.8. ")
    scree_analysis.add_argument("--unique", dest = "unique", action = "store_true", help = "Only retain cells with the top target that passed the freq_percent, which means each cell will be assigned with a unique sgRNA. ")
    scree_analysis.add_argument("--project", dest = "project", default = "perturb", type = str, help = "Project name. Default is 'perturb'. ")
    scree_analysis.add_argument("--species", dest = "species", default = "Hs", choices = ["Hs", "Mm"], help = "Species source of input data, one of 'Hs', 'Mm'. Defaule is 'Hs'. ")
    scree_analysis.add_argument("--fragments", dest = "fragments", required = False, help = "Directory to the fragments file. Only required for scATAC-seq based input. ")
    scree_analysis.add_argument("--prefix", dest = "prefix", default = ".", type = str, help = "Path to the results generated by SCREE. Default is current directory. ")
    scree_analysis.add_argument("--label", dest = "label", default = "", type = str, help = "The prefix label of previous output file. Notably, there needs a separator between default file names and the label, so label would be better to be like 'label_'. Default is ''. ")
#    scree_analysis.add_argument("--bar_width", dest = "bar_width", type = float, help = "Bar width of the barplot. Default is NULL, set to 90% of the resolution of the data. ")
    scree_analysis.add_argument("--feature_frac", dest = "feature_frac", default = 0.01, type = float, help = "A paramter for filtering low expressed genes or low accessible peaks. By default, only genes/peaks that have counts in at least that fractions of cells are kept. Default is 0.01. ")
    scree_analysis.add_argument("--nFeature", dest = "nFeature", default = 200, type = int, help = "Minimal count numbers in each cell. Default is 200.")
    scree_analysis.add_argument("--nCount", dest = "nCount", default = 1000, type = int, help = "Minimal detected feature numbers in each cell. Default is 1000. ")
    scree_analysis.add_argument("--percent.mt", dest = "percent.mt", default = 10, type = int, help = "Maximum mitochondrial gene percentage of each cell. This parameter can also be used for scATAC-seq based data to represent the minimum fraction of reads in peaks(FRiP). Default is 10(means 10%%). ")
    scree_analysis.add_argument("--blank", dest = "blank", action = "store_true", help = "Use blank control as negative control. With this parameter, cells assigned with no sgRNA will be removed. ")
    scree_analysis.add_argument("--NTC", dest = "NTC", default = "NTC", type = str, help = "The name of negative controls. Default is 'NTC'. ")
    scree_analysis.add_argument("--sg.split", dest = "sg.split", default = "_sgRNA", type = str, help = "String to split sgRNA name. Default is '_sgRNA', which means sgRNA named in the format like 'gene_sgRNA1'. ")
    scree_analysis.add_argument("--raster", dest = "raster", action = "store_true", help = "Convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. ")
    scree_analysis.add_argument("--ref_version", dest = "ref_version", default = "v75", type = str, choices = ["v75", "v86", "v79"], help = "Version of the reference genome(Ensembl), one of 'v75'(hg19/mm9), 'v79'(hg38/mm10), 'v86'(hg38). Default is 'v75'. ")
    scree_analysis.add_argument("--gene_type", dest = "gene_type", default = "Symbol", type = str, choices = ["Symbol", "Ensembl"], help = "Type of gene names, can be one of 'Symbol', 'Ensembl'. Default is 'Symbol'. ")
    scree_analysis.add_argument("--pval_cut", dest = "pval_cut", default = 0.05, type = float, help = "P-value cutoff of improved_scmageck_lr results. Default is 0.05. ")
    scree_analysis.add_argument("--score_cut", dest = "score_cut", default = 0.2, type = float, help = "Score cutoff of improved_scmageck_lr results. Default is 0.2. ")
    scree_analysis.add_argument("--article", dest = "article", default = "", type = str, help = "The link to the article of the dataset. ")
    scree_analysis.add_argument("--article_name", dest = "article_name", default = "", type = str, help = "Name of the article. ")
    scree_analysis.add_argument("--data", dest = "data", default = "", type = str, help = "Link to the data source. ")
    scree_analysis.add_argument("--data_name", dest = "data_name", default = "", type = str, help = "Name of the data, usually the GSE accession number. ")
    scree_analysis.add_argument("--padj_cut", dest = "padj_cut", default = 0.05, type = float, help = "Maximum adjust p_value. Default is 0.05. ")
    scree_analysis.add_argument("--logFC_cut", dest = "logFC_cut", default = 0.25, type = float, help = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25. ")
    
    args = parser.parse_args()
    
#     para.add_option("--datatype", default = "RNA", type = str, help = "Input data type, can be one of RNA/ATAC, only support scATAC-seq or scRNA-seq.")
#     para.add_option("--filetype", default = "multi", type = str, help = "Input file type, 'multi' indicates sgRNA sequence and mRNA/DNA sequence are not in the same files, 'single' indicates sgRNA sequence and mRNA sequence are in the same files. Notably, for scATAC-seq input, SCREE only support 'multi'.")
# #    para.add_option("--fasta", default = "", type = str, help = "File path of reference in FASTA format.")
# #    para.add_option("--gtf", "-g", default = "RNA", type = str, help = "Input data type, can be one of #RNA/ATAC, only support scATAC-seq or scRNA-seq.")
#     para.add_option("--reference", default = "", dest = "reference", type = str, help = "Reference path.")
#     para.add_option("--config", default = "", dest = "config", type = str, help = "Config file to generate reference index.")
#     para.add_option("--mkref", action = "store_true", dest = "mkref", help = "Generate reference index or not.")
#     para.add_option("--input", default = "", dest = "input", type = str, help = "File path of all the fastq files.")
#     para.add_option("--sample", default = "", dest = "sample", type = str, help = "Sample prefix.")
#     para.add_option("--output", default = ".", dest = "output", type = str, help = "Output folder to create.")
#     para.add_option("-m", dest = "memory", default = 50, type = int, help = "Maximum local memory to use. Unit: GiB")
#     para.add_option("-c", dest = "cores", default = 8, type = int, help = "Maximum local cores to use. ")
#     para.add_option("--csv", default = "", dest = "csv", type = str, help = "Config file to process 'multi' file type.")
#     para.add_option("--bin", action = "store_true", dest = "bin", help = "Generate a bin-based matrix instead of a peak-based matrix.")
#     para.add_option("--binsize", dest = "binsize", type = int, help = "Size of the genome bins to use. Only available with the parameter '--bin'.")
#     para.add_option("-n", "--process_n", dest = "process_n", type = int, help = "Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory. Only available with the parameter '--bin'.")
#     para.add_option("-l", "--chrLength", dest = "chrLength", type = str, help = "Path to the file which stores the chromosome length of the reference genome, always named as 'chrNameLength.txt'. Only available with the parameter '--bin'.")
#     para.add_option("--para", default = "", dest = "para", type = str, help = "Other parameters passed into cellranger or cellranger-atac, should be in string format like '--parameter1 value1 --parameter2 value2'")
    # para.add_option("--sgRNAgenomeGenerate", default = "", type = str, help = "File path of reference in FASTA format.")
    #options,args=para.parse_args()
    
    if args.version:
        print(version)
        exit(0)
    
    if args.subcommand == "preprocess":
        if len(sys.argv) == 2:
            scree_preprocess.print_help()
            exit(1)
        if args.filetype == "multi" and args.datatype == "RNA":
            script = "cellranger multi --id " + args.output + " --csv " + args.csv + " --localcores " + str(args.cores) + " --localmem " + str(args.memory) + " " + args.para
            os.system(script)
     #       exit(1)
        elif args.filetype == "single" and args.datatype == "RNA":
            script = "cellranger count --id " + args.output + " --transcriptome " + args.reference + " --fastqs " + args.input + " --sample " + args.sample + " --localcores " + str(args.cores) + " --localmem " + str(args.memory) + " " + args.para
            os.system(script)
     #       exit(1)
        elif args.datatype == "ATAC":
            if args.mkref:
                script = "cellranger-atac mkref --config " + args.config
                os.system(script)
            else:
                script = "cellranger-atac count --reference " + args.reference + " --fastqs " + args.input + " --sample " + args.sample + " --localcores " + str(args.cores) + " --localmem " + str(args.memory) + " --id " + args.output + " " + args.para
                os.system(script)
                if args.bin:
                    fragments = os.path.join(args.output, "outs", "fragments.tsv.gz")
                    out_file = os.path.join(args.output, "outs", "bin_matrix.rds")
                    Rscript_path = sys.path[0]
                    Rscript_path = os.path.join(Rscript_path, "bin.R")
                    os.system("Rscript " + Rscript_path  + " " + fragments + " " + str(args.binsize) + " " + str(args.process_n) + " " + args.chrLength + " " + out_file)
                                            
    elif args.subcommand == "analysis":
        if len(sys.argv) == 2:
            scree_analysis.print_help()
            exit(1)
        Rscript_path = sys.path[0]
        Rscript_path = os.path.join(Rscript_path, "analysis.R")
        script = "Rscript " + Rscript_path + " " + str(args.replicate) + " " + args.mtx + " " + args.sg_lib + " " + args.sg_format + " " + str(args.freq_cut) + " " + str(args.freq_percent) + " " + str(args.unique) + " " + args.project + " " + args.species + " " + args.fragments + " " + args.prefix + " " + args.label + " " + " " + str(args.feature_frac) + " " + str(args.nFeature) + " " + str(args.nCount) + " " + str(args.percent.mt) + " " + str(args.blank) + " " + args.NTC + " " + args.sg.split + " " + str(args.raster) + " " + args.ref_version + " " + args.gene_type + " " + str(args.pval_cut) + " " + str(args.score_cut) + " " + args.article + " " + args.article_name + " " + args.data + " " + args.data_name + " " + str(args.padj_cut) + " " + str(args.logFC_cut)
        os.system(script)
    else:
        parser.print_help()
        exit(1)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)