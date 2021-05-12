"""
Analysis pipeline for performing quality control, rRNA depletion, 
alignment, and feature counting for raw RNA-seq data for the manuscript, 
"Novel Classification of Human Ulcers using Host Gene Expression". 
Functionality is described in README.

Author: Blaine Gabriel Fritz 

Acknowledgements: Base code for this script was contributed by
                  Elio Rossi(Elio.Rossi@unimi.it)
"""


import argparse
import os
import sys
import subprocess
import shlex
import gzip 
import pdb


class RNAseqPipeline(object):
    def __init__(self, project_directory, data_directory, sample_name, output_directory, read1, 
                read2, reference_genome, annotation_gtf, rRNA_refs, threads, endedness):
        self.project_dir = project_directory
        self.sample_name = sample_name
        self.output_dir = output_directory
        self.data_directory = data_directory
        self.reference_genome = reference_genome
        self.annotation_gtf = annotation_gtf
        self.rRNA_refs = rRNA_refs
        self.threads = threads
        self.endedness = endedness

        self.read1 = read1
        self.read1_path = os.path.join(self.data_directory, self.read1)
        self.read1_trimmed = f"{self.sample_name}.trimmed_R1.fq.gz"
        self.read1_fq = f"{self.sample_name}.trimmed_R1.fq"
        self.rRNA_dep = f"{self.sample_name}.ribodep.fq"
        self.rRNA_reads = f"{self.sample_name}.rRNAreads.fq"
        self.sortmerna_R1 = f"{self.sample_name}_R1.ribodep.fastq.gz"
        self.bwa_aligned_bam = f"{self.sample_name}.aligned.bam"
        self.counts_outfile = f"{self.sample_name}.counts.txt"
        self.kraken_out = f"{self.sample_name}.kraken.output"
        self.kraken_report = f"{self.sample_name}.kraken.report"

        if self.endedness == "PE": 

            self.read2 = read2
            self.read2_path = os.path.join(self.data_directory, self.read2)
            self.read1_trimmed = f"{self.sample_name}.trimmed_R1.fq.gz" 
            self.read2_trimmed = f"{self.sample_name}.trimmed_R2.fq.gz"
            self.interleaved_reads = f"{self.sample_name}.interleaved.fq"
            self.rRNA_dep = f"{self.sample_name}.ribodep.intlvd.fq"
            self.rRNA_reads = f"{self.sample_name}.rRNAreads.intlvd.fq"
            self.sortmerna_R2 = f"{self.sample_name}_R2.ribodep.fastq.gz"       



    def trim_reads(self):
        #Define adapter sequences
        read1_adapt = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        read2_adapt = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

        #Run cutadapt
        if self.endedness == "PE":

            cutadapt_cmd = f"cutadapt -j 0 -a {read1_adapt} " \
                f"-A {read2_adapt} -m 20 " \
                f"-o {os.path.join(self.output_dir, self.read1_trimmed)} " \
                f"-p {os.path.join(self.output_dir, self.read2_trimmed)} " \
                f"{os.path.join(self.data_directory,self.read1_path)} " \
                f"{os.path.join(self.data_directory,self.read2_path)}"

        elif self.endedness == "SE":
            
            cutadapt_cmd = f"cutadapt -j 0 -a {read1_adapt} -m 20 " \
                f"-o {os.path.join(self.output_dir, self.read1_trimmed)} " \
                f"{os.path.join(self.data_directory,self.read1_path)} " 

        trm = subprocess.Popen(
            shlex.split(cutadapt_cmd), 
            shell=False,
            stdout=open(
                f"{os.path.join(self.output_dir,self.sample_name)}_cutadapt.log",'w')
            )  

        trm.communicate()

        if trm.returncode != 0:
            print("ERROR: Trimming Failed")
            

    def deplete_rRNA(self):

        if self.endedness == "PE":
            #Interleave the reads from the two, trimmed read files
            with open(os.path.join(self.output_dir,self.interleaved_reads), "w") as intlvd_fq:

                intlvd_cmd = f"seqtk mergepe {os.path.join(self.output_dir, self.read1_trimmed)} " \
                    f"{os.path.join(self.output_dir, self.read2_trimmed)} "

                seqtk = subprocess.Popen(shlex.split(intlvd_cmd), shell=False,
                                         stdout=intlvd_fq,
                                         )

                seqtk.communicate()

                if seqtk.returncode != 0:
                    print("Error while interleaving reads. Exiting...")

            #SortMeRNA Command

            rrna_cmd = "sortmerna " \
                f"--ref {self.rRNA_refs} " \
                f"--reads {os.path.join(self.output_dir, self.interleaved_reads)} " \
                f"--paired_in " \
                f"-a {self.threads} " \
                f"--log --fastx " \
                f"--aligned {os.path.join(self.output_dir,self.rRNA_reads).replace('.fq', '')} " \
                f"--other {os.path.join(self.output_dir,self.rRNA_dep).replace('.fq', '')} " \
                f"-v"
            
            srtmrna = subprocess.Popen(shlex.split(rrna_cmd), shell=False)
            
            srtmrna.communicate()

            if srtmrna.returncode != 0:
                    print("ERROR while removing rRNA")

            #R1
            seqtk_cmd1 = f"seqtk seq -1 " \
                f"{os.path.join(self.output_dir, self.rRNA_dep)}"
            
            #R2    
            seqtk_cmd2 = f"seqtk seq -2 " \
                f"{os.path.join(self.output_dir, self.rRNA_dep)}"

            gzip_cmd = "gzip"
            
            
            with open(os.path.join(self.output_dir, self.sortmerna_R1), "wb") as r1_f:    
                seqtk_1 = subprocess.Popen(
                    shlex.split(seqtk_cmd1), shell = False, 
                    stdout=subprocess.PIPE)
                
                gzip = subprocess.Popen(
                        shlex.split(gzip_cmd),
                        stdin=seqtk_1.stdout,
                        stdout=r1_f
                    )

                gzip.communicate()
                
                if gzip.returncode != 0:
                    print("Error: Something went wrong with rRNA depletion")

            with open(os.path.join(self.output_dir, self.sortmerna_R2), "wb") as r2_f:
                
                seqtk_2 = subprocess.Popen(
                    shlex.split(seqtk_cmd2), shell = False, 
                    stdout=subprocess.PIPE)
                
                gzip = subprocess.Popen(
                        shlex.split(gzip_cmd),
                        stdin=seqtk_2.stdout,
                        stdout=r2_f
                    )

                gzip.communicate()
                
                if gzip.returncode != 0:
                    print("Error: Something went wrong with rRNA depletion")              


        if self.endedness == "SE":
            
            unzip_cmd = f"gunzip {os.path.join(self.output_dir, self.read1_trimmed)}"

            unzip = subprocess.Popen(shlex.split(unzip_cmd), shell=False)

            unzip.communicate()

            if unzip.returncode != 0:
                print("ERROR unzipping trimmed reads")

            rrna_cmd = "sortmerna " \
                f"--ref {self.rRNA_refs} " \
                f"--reads {os.path.join(self.output_dir, self.read1_fq)} " \
                f"-a {self.threads} " \
                f"--log --fastx " \
                f"--aligned {os.path.join(self.output_dir,self.rRNA_reads).replace('.fq', '')} " \
                f"--other {os.path.join(self.output_dir,self.rRNA_dep).replace('.fq', '')} " \
                f"-v"
            
            srtmrna = subprocess.Popen(shlex.split(rrna_cmd), shell=False)

            srtmrna.communicate()
            
            if srtmrna.returncode != 0:
                print("ERROR: SortmeRNA FAILED. ")
            

            with open(os.path.join(self.output_dir, self.sortmerna_R1), "wb") as r1_f:
            
                gzip_cmd = f"gzip -c {os.path.join(self.output_dir, self.rRNA_dep)}"

                gzip = subprocess.Popen(
                        shlex.split(gzip_cmd),
                        stdout=r1_f, 
                        shell = False
                    )

                gzip.communicate()

                if gzip.returncode != 0:
                    print("ERROR zipping sortmeRNA output")          

                      
    def run_alignment(self):

        #Define the commands 
        if self.endedness =="PE":

            bwa_cmd = f"bwa mem -t {self.threads} " \
                f"{self.reference_genome} " \
                f"{os.path.join(self.output_dir, self.sortmerna_R1)} " \
                f"{os.path.join(self.output_dir, self.sortmerna_R2)}"

        
        if self.endedness == "SE":

            bwa_cmd = f"bwa mem -t {self.threads} " \
                f"{self.reference_genome} " \
                f"{os.path.join(self.output_dir, self.sortmerna_R1)} " 

       
        smt_cmd = f"samtools view -Sbh -@ {self.threads}"    
        
        with open(os.path.join(self.output_dir, self.bwa_aligned_bam), "wb") as algn_file:

            bwa = subprocess.Popen(shlex.split(bwa_cmd), shell=False,
                                   stdout=subprocess.PIPE)

            smt = subprocess.Popen(shlex.split(smt_cmd), stdin=bwa.stdout,
                                   stdout=algn_file)

            smt.communicate()

            if smt.returncode != 0:
                print("ERROR: SortmeRNA FAILED. ")
                    

    def count_reads(self):

        if self.endedness == "PE":
            count_cmd = f"featureCounts -p -T {self.threads} -t exon -g gene " \
                f"-o {os.path.join(self.output_dir, self.counts_outfile)} " \
                f"-a {self.annotation_gtf} " \
                f"{os.path.join(self.output_dir, self.bwa_aligned_bam)}"

        if self.endedness == "SE":
            count_cmd = f"featureCounts -T {self.threads} -t exon -g gene " \
                f"-o {os.path.join(self.output_dir, self.counts_outfile)} " \
                f"-a {self.annotation_gtf} " \
                f"{os.path.join(self.output_dir, self.bwa_aligned_bam)}"
        

        count = subprocess.Popen(shlex.split(count_cmd), shell = False)

        count.communicate()

        if count.returncode != 0:
            print("ERROR: featureCounts FAILED. ")
            

    def zip_files(self): 

        if self.endedness == "PE": 

            pigz_cmd = f"pigz {os.path.join(self.output_dir, self.interleaved_reads)} " \
                f"{os.path.join(self.output_dir, self.rRNA_dep)} " \
                f"{os.path.join(self.output_dir, self.rRNA_reads)}"
        
        if self.endedness == "SE":
            pigz_cmd = f"pigz {os.path.join(self.output_dir, self.read1_fq)} " \
                f"{os.path.join(self.output_dir, self.rRNA_dep)} " \
                f"{os.path.join(self.output_dir, self.rRNA_reads)}" 

            
        pigz = subprocess.Popen(shlex.split(pigz_cmd), shell = False)
        
        pigz.communicate()

        if pigz.returncode != 0:
         print("ERROR: Pigz failed to zip the files properly. ")

    def run_kraken2(self): 

        if self.endedness == "PE":

            cmd_str = f"kraken2 --db /your/path/here/data/databases/kraken2/kraken2_std_db --gzip-compressed " \
                f"--threads {self.threads} " \
                f"--output - " \
                f"--report {os.path.join(self.output_dir, self.kraken_report)} " \
                f"--paired {self.read1_path}  {self.read2_path}"

        if self.endedness == "SE":

            cmd_str = f"kraken2 --db /your/path/here/data/databases/kraken2/kraken2_std_db --gzip-compressed " \
                f"--threads {self.threads} " \
                f"--output - " \
                f"--report {os.path.join(self.output_dir, self.kraken_report)} " \
                f"{self.read1_path}"


        cmd = subprocess.Popen(shlex.split(cmd_str),
                   shell=False)
        
        cmd.communicate()

        if cmd.returncode != 0:
            print("ERROR: KRAKEN2 didn't work properly.")


def main():
    """ 1. Trim reads
        2. Remove rRNA
        3. Align to reference genome
        4. Count reads
    """

    #Parse input args and create a dictionary
    cmd_args = parse_args()
    cmd_args_dict = vars(cmd_args)

    # Change to the project directory

    analysis = RNAseqPipeline(**cmd_args_dict)
    
    analysis.trim_reads()
    analysis.deplete_rRNA()
    analysis.run_alignment()
    analysis.count_reads()
    analysis.zip_files()
    analysis.run_kraken2()



def parse_args():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("-p", "--project_directory",
                            help="Path to main project directory",
                            required=True
                            )
    arg_parser.add_argument("-e", "--endedness", 
                            choices=['SE','PE'], 
                            help="Paired-end or single-end reads",
                            required=True
                            )
    arg_parser.add_argument("-d", "--data_directory",
                            help="Path to data directory",
                            required=True
                            )
    arg_parser.add_argument("-s", "--sample_name",
                            help="Name of Sample. Taken from metadata",
                            required=True
                            )                        
    arg_parser.add_argument("-o", "--output_directory",
                            help="Path to output folder",
                            required=True
                            )
    arg_parser.add_argument("-r1", "--read1",
                            help="Name of forward read. Taken from metadata",
                            required=True
                            )
    arg_parser.add_argument("-r2", "--read2",
                            help="Name of forward read (if PE). Taken from metadata",
                            required=False
                            )
    arg_parser.add_argument("-ref", "--reference_genome",
                            help="Path to main project directory",
                            required=True
                            )
    arg_parser.add_argument("-a", "--annotation_gtf",
                            help="Path to annotation GTF file",
                            required=True
                            )
    arg_parser.add_argument("--rRNA_refs",
                        help="Paths to rrna depletion databases",
                        required=True
                        )
    arg_parser.add_argument("-t", "--threads",
                        help="Number of threads to use",
                        required=True
                        )
    return arg_parser.parse_args()


if __name__ == '__main__':
    main()            