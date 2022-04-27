from Bio import SeqIO
import argparse
import os
import subprocess
import io
import multiprocessing as mp
import shutil
import pandas as pd
import glob

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def create_folder(mypath):

    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return

###########################################################

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--in_fasta_file",
                        help="input_fasta", 
                        type=str, 
                        default=snakemake.input.fasta)
    parser.add_argument("--in_tab_file", 
                        help="input_table", 
                        type=str, 
                        default=snakemake.params.viral_affi)
    parser.add_argument("--min_contig_size", 
                        help="minimum contig size", 
                        type=int, 
                        default=snakemake.params.cutoff)
    parser.add_argument('--outfile', 
                        help='output_folder', 
                        type=str, 
                        default=snakemake.output.tsv)
    parser.add_argument('--out_faa', 
                        help='output_folder', 
                        type=str, 
                        default=snakemake.output.faa)
    parser.add_argument('--tmp_dir', 
                        help='tmp folder', 
                        type=str, 
                        default=os.path.join(snakemake.params.output_dir, "tmp"))
    parser.add_argument('--contigs_per_file', 
                        help='The number of contigs in each slice', 
                        type=int, 
                        default=1)
    parser.add_argument('--job_number', 
                        type=int, 
                        default=snakemake.threads)
    parser.add_argument('--threads_per_process', 
                        type=int, 
                        default=16)    
    # Parse arguments
    args = parser.parse_args()

    return args

###########################################################

def main(args):
    
    #split the input file
    record_iter = SeqIO.parse(open(args.in_fasta_file), "fasta")
    files_to_run = []


    for i, batch in enumerate(batch_iterator(record_iter, args.contigs_per_file)):
        try:
            group_name=f'group_{i+1}'
            dir_name=f'{args.tmp_dir}/{group_name}'
            filename = f"{dir_name}/{group_name}.fasta"

            create_folder(f"{args.tmp_dir}/{group_name}")

            with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
                print(f"Wrote {count} records to {filename}")
                files_to_run.append((group_name, dir_name, filename))
        except IOError as ioerror:
            print(ioerror)

    pool = mp.Pool(args.job_number)
    results = pool.map(run_job, files_to_run)
    pool.close()

    df = pd.concat(results)
    df.to_csv(args.outfile, sep='\t')

    merge_fasta(args.tmp_dir, args.out_faa)

    print(f'{bcolors.OKBLUE} ------ DONE! ----------- {bcolors.ENDC}')
    # shutil.rmtree(args.tmp_dir)

###########################################################

def merge_fasta(tmp_dir, out_fasta):
    all_prot = glob.glob(tmp_dir, "*", "dram-v-output", "*.faa")

    with open(out_fasta, "wt") as w_file:
        for faa_file in all_prot:
            parser = SeqIO.parse(faa_file, "fasta")
            for prot in parser:
                prot.name = prot.description = ''
                SeqIO.write(prot, w_file, "fasta")

    return

###########################################################

def run_job(group_tuple):
    if args.in_tab_file:
        job_str=f'DRAM-v.py annotate -i {group_tuple[2]} -v {args.in_tab_file} -o {group_tuple[1]}/dram-v-output --skip_trnascan --threads {args.threads_per_process} --min_contig_size {args.min_contig_size}'
    else:
        job_str=f'DRAM-v.py annotate -i {group_tuple[2]} -o {group_tuple[1]}/dram-v-output --skip_trnascan --threads {args.threads_per_process} --min_contig_size {args.min_contig_size}'
    stdout, stderr = execute(job_str)
    print(f"----DRAMv - stdout----\n{stdout.decode('utf8')}\n----DRAMv - stderr----\n{stderr.decode('utf8')}\n")
    if os.path.isfile(f'{group_tuple[1]}/dram-v-output/annotations.tsv'):
        df = pd.read_csv(f'{group_tuple[1]}/dram-v-output/annotations.tsv', sep='\t', index_col=0)
    else:
        df = pd.DataFrame()
    return df

###########################################################

def execute(command):
    print(f'Executing {command}')
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    stdout, stderr = process.communicate()
    return stdout, stderr

###########################################################

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

###########################################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

###########################################################

if __name__ == "__main__":
    args = parse_args()
    main(args)
