import os
def main():
    from optparse import OptionParser
    para = OptionParser(usage = 'SCREE [--help/-h][--datatype][--filetype][--reference][--config][--mkref][--input][--sample][--output][-m][-c][--csv]')
    para.add_option("--datatype", default = "RNA", type = "string", help = "Input data type, can be one of RNA/ATAC, only support scATAC-seq or scRNA-seq.")
    para.add_option("--filetype", default = "multi", type = "string", help = "Input file type, 'multi' indicates sgRNA sequence and mRNA/DNA sequence are not in the same files, 'single' indicates sgRNA sequence and mRNA sequence are in the same files. Notably, for scATAC-seq input, SCREE only support 'multi'.")
#    para.add_option("--fasta", default = "", type = "string", help = "File path of reference in FASTA format.")
#    para.add_option("--gtf", "-g", default = "RNA", type = "string", help = "Input data type, can be one of #RNA/ATAC, only support scATAC-seq or scRNA-seq.")
    para.add_option("--reference", default = "", type = "string", help = "Reference path.")
    para.add_option("--config", default = "", type = "string", help = "Config file to generate reference index.")
    para.add_option("--mkref", action = "store_true", help = "Generate reference index or not.")
    para.add_option("--input", default = "", type = "string", help = "File path of all the fastq files.")
    para.add_option("--sample", default = "", type = "string", help = "Sample prefix.")
    para.add_option("--output", default = "./", type = "string", help = "Output folder to create.")
    para.add_option("-m", dest = "memory", default = 50, type = "int", help = "Maximum local memory to use. Unit: GiB")
    para.add_option("-c", dest = "cores", default = 8, type = "int", help = "Maximum local cores to use. ")
    para.add_option("--csv", default = "", type = "string", help = "Config file to process 'multi' file type.")
    para.add_option("--para", default = "", type = "string", help = "Other parameters passed into cellranger or cellranger-atac, should be in string format like '--parameter1 value1 --parameter2 value2'")
    # para.add_option("--sgRNAgenomeGenerate", default = "", type = "string", help = "File path of reference in FASTA format.")
    options,args=para.parse_args()
    
    if options.filetype == "multi" and options.datatype == "RNA":
        script = "cellranger multi --id " + options.output + " --csv " + options.csv + " --localcores " + options.cores + " --localmem " + options.memory + " " + options.para
        os.system(script)
        exit(1)
    elif options.filetype == "single" and options.datatype == "RNA":
        script = "cellranger count --id " + options.output + " --transcriptome " + options.reference + " --fastqs " + options.input + " --sample " + options.sample + " --localcores " + options.cores + " --localmem " + options.memory + " " + options.para
        os.system(script)
        exit(1)
    elif options.datatype == "ATAC":
        if options.mkref:
            script = "cellranger-atac mkref --config " + options.config
            os.system(script)
        else:
            script = "cellranger-atac count —-reference " + options.reference + " --fastqs " + options.input + " --sample " + options.sample + " --localcores " + options.cores + " --localmem " + options.memory + " --id " + options.output + " " + options.para
            os.system(script)
        exit(1)
    
    exit(0)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)