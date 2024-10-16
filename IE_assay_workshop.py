#!/usr/bin/env python3

# Copyright 2024 Dorian J. Feistel
# Author: Dorian J. Feistel
# Email: djfeistel@gmail.com
# Date: 2024.06.12
# License: GPLv3 <http://www.gnu.org/licenses/gpl-3.0.html>
# Description: Script for in silico Inclusivity/Exclusivity Assay

'''
things to include:
    coverage plot: qould this be for all hits or for just the max coverage per query ?
    
'''
import os
import sys
import subprocess
import logging
import shutil
import argparse
import pandas as pd
import numpy as np
import tempfile
import random, string
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from statsmodels.stats.proportion import proportions_ztest


def setup_logging_info(logname, verbose=False):
    logname = "IE_assay.log" if logname is None else logname + ".log"
    if os.path.isfile(logname):
        os.remove(logname)
    # Create the handlers list with FileHandler
    handlers = [logging.FileHandler(logname)]
    
    if verbose:
        handlers.append(logging.StreamHandler())

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )
    return logname

def logging_args(args):
    vargs = vars(args)
    for arg, val in vargs.items():
        logging.info(f"Parameters: {arg} = {val}")

def check_file_existance(infile:str):
    if not os.path.isfile(infile):
        logging.error(f"Input file '{infile}' not found.")
        raise FileNotFoundError(f"Input file '{infile}' not found.")
    logging.info(f"File existance check passed: {infile} is a real file")

def check_fasta_existance(infasta: str):
    '''check all records in a fasta file more efficiently'''
    try:
        # Try parsing the FASTA file using SeqIO
        records = list(SeqIO.parse(infasta, "fasta"))
        # If parsing succeeds and records are found, the file is correctly formatted
        if records:
            logging.info(f"Fasta check passed: {len(records)} records found.")
        else:
            logging.warning(f"No records found in the file: {infasta}.")
            raise ValueError(f"No records found in the file: {infasta}.")
    except Exception as e:
        logging.error(f"Error parsing FASTA file '{infasta}': {str(e)}")
        raise ValueError(f"Error parsing FASTA file '{infasta}': {str(e)}")
                             
def is_blastn_available():
    try:
        # Run 'blastn' with the '--version' flag
        result = subprocess.run(['blastn', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Check if the command was successful
        if result.returncode != 0:
            logging.error("blastn is not available.")
            raise RuntimeError("blastn is not available.")
    except FileNotFoundError:
        logging.error("blastn command not found.")
        raise RuntimeError("blastn command not found.")

def load_primer_assays(primer_file:str): 
    names = ['primer_name', 'primer_orientation', 'primer_assay', 'sequence']
    df = pd.read_csv(primer_file, sep='\t', header=None, names=names)

    combinations = []
    # primer_dict = df.groupby('primer_assay')['primer_name'].apply(list).to_dict()
    # print(primer_dict)
    grouped = df.groupby('primer_assay')
    primer_dict = {}
    for name, group in grouped:

        forward = group[group['primer_orientation'] == '+']['primer_name']
        reverse = group[group['primer_orientation'] == '-']['primer_name']
        probe = group[group['primer_orientation'] == '.']['primer_name']
        
        primer_dict[name] = [forward, reverse, probe]
        comb = list(product(forward, reverse, probe))
        combinations.extend(comb)
    print(primer_dict)
    return combinations

class DegenerateSequenceProcessor:
    def __init__(self, bedfile, dir):
        self.bed = bedfile #infile
        self.bed_df = None # pandas df
        self.dir = dir # temp dir
        self.sequences = None # biopython seqs
        self.oligo_dict = {}
        self.combinations = []
        self.seq_records = {}
        self.output_fasta = "Disambiguated_sequences.fasta"
        self.degenerate_bases = {
            'R': 'AG',  # Purine
            'Y': 'CT',  # Pyrimidine
            'S': 'GC',  # Strong Interaction
            'W': 'AT',  # Weak Interaction
            'K': 'GT',  # Keto
            'M': 'AC',  # Amino
            'B': 'CGT', # Not A
            'D': 'AGT', # Not C
            'H': 'ACT', # Not G
            'V': 'ACG', # Not T
            'N': 'ACGT' # Any Nucleotide
        }
        
    def bed2df(self):
        names = 'primer_name primer_orientation primer_assay sequence'.split(" ")
        self.bed_df = pd.read_csv(self.bed, sep='\t', header=None, names=names)
        self.bed_df['primer_orientation'] = pd.Categorical(self.bed_df['primer_orientation'], categories=['+', '.', '-'], ordered=True)
        self.bed_df = self.bed_df.sort_values(by=['primer_assay', 'primer_orientation'])

    def primer_name_orientation_dict(self):
        return dict(zip(self.bed_df['primer_name'], self.bed_df['primer_orientation']))
    
    def primer_name_assay_dict(self):
        return dict(zip(self.bed_df['primer_name'], self.bed_df['primer_assay']))
    
    def primer_assay_name_dict(self):
        return  self.bed_df.groupby('primer_assay')['primer_name'].apply(list).to_dict()

    def bed2seqs(self):
        '''convert bed-like tsv file to fasta'''
        try:
            orient_description = {"+": "Forward", "-": "Reverse", ".": "Probe"}
            logging.info(f"Opening {self.bed}.")
            self.sequences = [
                SeqRecord(Seq(seq), 
                          id=primer_name, 
                          description=f"{primer_assay}###{orient_description[orientation]}###"
                          )
                    for primer_name, seq, primer_assay, orientation in zip(
                        self.bed_df['primer_name'], 
                        self.bed_df['sequence'], 
                        self.bed_df['primer_assay'], 
                        self.bed_df['primer_orientation']
                    )
            ]

            if self.sequences:
                logging.info(f"{os.path.basename(self.bed)} has a total of {len(self.sequences)} sequences.")
                
            else:
                logging.error(f"No sequences found in {self.bed}.")
                raise ValueError(f"No sequences found in {self.bed}.")

        except FileNotFoundError as fnf_error:
            logging.error(fnf_error)
        except ValueError as val_error:
            logging.error(val_error)
        except Exception as e:
            logging.error(f"An unexpected error occurred: {e}")

    def generate_sequence_combinations(self, seq):
        '''
        Generates all possible sequences by replacing each degenerative base 
        in the input sequence with its corresponding nucleotides
        '''
        try:
            # Generate all possible replacements for each base in the sequence
            possible_replacements = [self.degenerate_bases.get(base, base) for base in seq]
            # Generate all possible combinations
            self.combinations = [''.join(p) for p in product(*possible_replacements)]
        
        except Exception as e:
            logging.error(f"An error occurred while generating sequence combinations: {e}")

    def process_degenerate_sequences(self):
        
        for seq_record in self.sequences:
            
            id = seq_record.id.strip()
            seq = str(seq_record.seq).upper().strip()
            primer_assay = seq_record.description.strip().split("###")[0]
           
            orientation = seq_record.description.strip().split("###")[1]

            if primer_assay not in self.oligo_dict:
                self.oligo_dict[primer_assay] = {}
            if orientation not in self.oligo_dict[primer_assay]:
                self.oligo_dict[primer_assay][orientation] = []

            # Check for invalid characters
            if any(base not in "AGCTRYWSKMBDHVN" for base in seq):
                logging.error(f"Sequence {id} contains invalid base characters.")
                raise ValueError(f"Sequence {id} contains base invalid characters.")
            
            try:
                self.generate_sequence_combinations(seq)

            except Exception as e:
                logging.error(f"An error occurred while generating combinations for sequence {id}: {e}")
                raise ValueError(f"An error occurred while generating combinations for sequence {id}: {e}")
    
            for idx, comb in enumerate(self.combinations, start=1):

                new_header = f"{id}_oligo{idx}"
                new_seq_record = SeqRecord(Seq(comb), id=new_header, description=f"Oligo {idx} out of {len(self.combinations)}; part of primer assay {primer_assay}; Orientaiton is {orientation}")
                self.seq_records[new_header] = new_seq_record
                logging.info(f"Added new sequence record: {new_header}")
                
                # used later in blaststatistics
                self.oligo_dict[primer_assay][orientation].append(new_header)

    def write_degenerate_sequences(self):
        output_pathway = os.path.join(self.dir, self.output_fasta)

        try:

            SeqIO.write(list(self.seq_records.values()), output_pathway, "fasta")
            logging.info(f"Generated {len(self.seq_records)} sequences and saved to {self.output_fasta}")
        except Exception as e:
            logging.error(f"Error writing output file {self.output_fasta}: {e}")

class BlastUtilities:
    def __init__(
            self, 
            query: str, 
            reference: str, 
            dir: str, 
            assay: str, 
            output: str,
            threads: str = "1"
        ):

        self.query = query
        self.reference = reference
        self.dir = dir
        self.assay = assay
        self.output = output
        self.threads = threads
        self.db = None
        self.blast_file = None

    def is_blastn_available(self):
        try:
            # Run 'blastn' with the '--version' flag
            result = subprocess.run(['blastn', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Check if the command was successful
            if result.returncode != 0:
                logging.error("blastn is not available.")
                raise RuntimeError("blastn is not available.")
        except FileNotFoundError:
            logging.error("blastn command not found.")
            raise RuntimeError("blastn command not found.")
        logging.info(f"blastn test suscessful.")

    def extract_reference_deflines(self):
        '''extract all reference sequence IDs (sseqid) from referecne fasta file'''
        deflines = []
        for record in SeqIO.parse(self.reference, "fasta"):
            deflines.append(record.id)
        return deflines

    def make_blast_db(self) -> None:
        '''Create BLAST database in nucleotide format'''

        # check if db already exists, if so do not create it but use it
        db_extensions =  ['.nhr', '.nin', '.nsq']
        for file in os.listdir(os.path.dirname(os.path.abspath(self.reference))):
            
            if file.startswith(os.path.basename(self.reference)) and any(file.endswith(ext) for ext in db_extensions):
                self.db = self.reference
                logging.info(f"BLAST databse for {os.path.basename(self.reference)} already exists in {os.path.dirname(self.reference)}.")
                logging.info(f"No BLAST db will be included in the final output. Using db={os.path.abspath(self.db)}")
                return
           
        blast_db = os.path.join(self.dir, "db")
        self.db = os.path.join(blast_db, os.path.basename(self.reference))
        command = ['makeblastdb', '-in', self.reference, '-dbtype', 'nucl', '-out', self.db]

        try:
            logging.info(f"Running makblastdb command: {' '.join(command)}. Additional time may be needed to create DB for larger files.")
            result = subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.stdout:
                logging.info(result.stdout)
            logging.info("makeblastdb completed successfully.")
        except subprocess.CalledProcessError as e:
            logging.error(f"makeblastdb failed with return code {e.returncode}.")
            raise RuntimeError(f"makeblastdb failed with return code {e.returncode}.")
        except FileNotFoundError:
            logging.error("makeblastdb command not found.")
            raise RuntimeError("makeblastdb command not found.")
        shutil.copy2(self.reference, blast_db)

    def run_blastn(self) -> str:
        '''Run BLASTn to compare query sequences against a BLAST database'''
        extention = "_raw_blast_results.tsv"
        output = f"{self.assay}{extention}" if self.output is None else f"{self.output}{extention}"

        self.blast_file = os.path.join(self.dir, output)
        outfmt = "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen qseq sseq sstrand"
        strict_parms = ['-word_size', '7', '-reward', '1', '-penalty', '-3', '-gapopen', '5', '-gapextend', '2', '-evalue', '1000']
        blast_command = ['blastn', '-num_threads', self.threads, '-query', os.path.join(self.dir, self.query), '-db', self.db, '-out', self.blast_file, '-outfmt', outfmt, '-max_target_seqs', '10000000']
        command  = blast_command + strict_parms
        try:
            logging.info(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            logging.info(f"BLASTn completed successfully. Output saved to {output}.")
            return output
        except subprocess.CalledProcessError as e:
            logging.error(f"BLASTn failed with return code {e.returncode}. Command: {' '.join(command)}")
            raise RuntimeError(f"BLASTn failed with return code {e.returncode}.")
        except FileNotFoundError:
            logging.error("BLASTn command not found.")
            raise RuntimeError("BLASTn command not found.")
        except Exception as e:
            logging.error(f"An unexpected error occurred: {str(e)}")
            raise RuntimeError(f"An unexpected error occurred: {str(e)}")
        
class BlastFilter:
    def __init__(
            self, 
            blast_file: str,
            output: str,
            dir: str, 
            assay: str,
            primer_name_orientation_dict: dict,
            primer_name_assay_dict: dict,
            seqid: float = None,
            coverage: float = None,
            mismatch: int = None,
            uncovered: int = None,
            snps: int = None,
        ):

        self.blast_file = blast_file # raw blast file
        self.raw_blast_data = pd.DataFrame()
        self.filtered_blast = pd.DataFrame()
        self.output = output
        self.assay = assay
        self.dir = dir
        self.primer_name_orientation_dict = primer_name_orientation_dict
        self.primer_name_assay_dict = primer_name_assay_dict
        self.seqid = seqid
        self.coverage = coverage
        self.mismatch = mismatch
        self.uncovered = uncovered
        self.snps = snps
        

    def load_assay_parameters(self):
        '''Set parameters when flags are i or e'''
        if self.assay == 'i':
            self.seqid = self.seqid if self.seqid is not None else None
            self.coverage = self.coverage if self.coverage is not None else None
            self.uncovered = self.uncovered if self.uncovered is not None else 2
            self.mismatch = self.mismatch if self.mismatch is not None else None
            self.snps = self.snps if self.snps is not None else 2
            # might want to add in addtional output data about same primer aligning to different locations

        elif self.assay == 'e':
            self.seqid = self.seqid if self.seqid is not None else 80
            self.coverage = self.coverage if self.coverage is not None else 0.8
            self.uncovered = self.uncovered if self.uncovered is not None else None
            self.mismatch = self.mismatch if self.mismatch is not None else None
            self.snps = self.snps if self.snps is not None else None

        params_to_log = {k: v for k, v in vars(self).items() if k in ['assay', 'seqid', 'coverage', 'uncovered', 'mismatch', 'snps']}
        logging.info(f"Assay parameters use to filter raw BLAST results: {params_to_log}")

    def load_data(self):
        '''Load raw BLAST data from its relative pathway.'''
        try:
            logging.info(f"Loading BLAST data from {os.path.basename(self.blast_file)}")

            columns = [
                "qseqid", "sseqid", "pident", "length", "mismatch",
                "gapopen", "gaps", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore", "qlen", "qseq", "sseq", "sstrand"
            ]
            self.raw_blast_data = pd.read_csv(self.blast_file, sep='\t', header=None, names=columns)
            self.raw_blast_data['coverage'] = round(self.raw_blast_data['length'] / self.raw_blast_data['qlen'], 3) # breadth of sequencing coverage
            self.raw_blast_data[['z_score', 'p_value']] = self.raw_blast_data.apply(self.proportion_test_alignment, snps=self.snps, axis=1)
            self.raw_blast_data.to_csv(self.blast_file, sep='\t', index=False)
            self.filtered_blast = self.raw_blast_data.copy()
            
            logging.info("BLAST data loaded successfully and sequencing breadth of coverage calculated (length/qlen).")

        except Exception as e:
            logging.error(f"Failed to load BLAST data: {str(e)}")
            raise RuntimeError(f"Failed to load BLAST data: {str(e)}")

    def filter_by_seqid(self):
        self.filtered_blast = self.filtered_blast[self.filtered_blast['pident'] >= self.seqid]
        logging.info(f"Filtered data by minimum percent identity of {self.seqid}% with {len(self.filtered_blast)} records remain.")
    
    def filter_by_coverage_breadth(self):
        self.filtered_blast = self.filtered_blast[self.filtered_blast['coverage'] >= self.coverage]
        logging.info(f"Filtered data by minimum coverage (length/qlen) of {self.coverage} with {len(self.filtered_blast)} records remain.")
        
    def filter_by_coverage_uncovered(self):
        self.filtered_blast = self.filtered_blast[(self.filtered_blast['qlen'] - self.filtered_blast['length']) <= self.uncovered]
        logging.info(f"Filtered data by maximum number of uncovered bases {self.uncovered}. {len(self.filtered_blast)} records remain.")

    def filter_by_mismatch(self):
        self.filtered_blast = self.filtered_blast[self.filtered_blastf['mismatch'] <= self.mismatch]
        logging.info(f"Filtered data by maximum mismatches of {self.mismatch}. {len(self.filtered_blast)} records remain.")

    def filter_by_SNPs(self):
        self.filtered_blast = self.filtered_blast[(self.filtered_blast['mismatch'] + self.filtered_blast['gaps']) <= self.snps]
        logging.info(f"Filtered data by maximum SNPs of {self.snps}. {len(self.filtered_blast)} records remain.")

    def filter_highest_bitscore_per_qseqid(self):
        max_bitscores = self.filtered_blast.groupby(['qseqid', 'sseqid'])['bitscore'].idxmax()
        self.filtered_blast = self.filtered_blast.loc[max_bitscores]
        logging.info(f"Filtered data to keep the best matching qseqid vs sseqid based on the highest bitscore. {len(self.filtered_blast)} records remain.")

    def alignment_distribution(self):
        logging.info("Creating alignment distribution on raw BLAST results.")
        logging.getLogger('matplotlib').setLevel(logging.WARNING)

        sns.set_theme(style="whitegrid")

        plt.figure(figsize=(10, 6))

        ax = sns.histplot(
            data=self.raw_blast_data, 
            x='length', 
            hue='qseqid', 
            element='step', 
            stat='frequency',
            multiple='layer',
            common_norm=False)

        ax.set_yscale('log')
        current_min, current_max = ax.get_ylim()
        ax.set_ylim(current_min, 10**((np.log10(current_max) + 1)))
        plt.xlabel('Alignment Length (bp)')
        plt.ylabel('Frequency (log scale)')
        plt.title('Frequency of Alignment Lengths by Oligos')

        output = f"{self.assay}_alignment_distribution.png" if self.output is None else f"{self.output}_alignment_distribution.png"
        output = os.path.join(self.dir, output)
        plt.savefig(output, dpi=150)

        logging.info(f"Alignment distribution plot saved to {output}")

    def alignment_boxplot(self):
        logging.getLogger('matplotlib').setLevel(logging.WARNING)

        # Group the data by 'qseqid' and get the 'length' for each group
        unique_qseqid = self.raw_blast_data['qseqid'].unique()
        grouped_data = [self.raw_blast_data[self.raw_blast_data['qseqid'] == qseqid]['length'] for qseqid in unique_qseqid]
        
        # Set up the figure
        plt.figure(figsize=(20, 6))
        
        # Create the box plot
        plt.boxplot(grouped_data, labels=unique_qseqid, patch_artist=True)
        
        # Overlay scatter plots for each qseqid
        for i, qseqid in enumerate(unique_qseqid):
            # Get the 'length' values for the current 'qseqid'
            data_points = self.raw_blast_data[self.raw_blast_data['qseqid'] == qseqid]['length']

            # Scatter plot of the actual data points with slight random jitter to prevent overlap
            jittered_x = np.random.normal(i + 1, 0.04, size=len(data_points))  # Add jitter to x positions
            plt.scatter(jittered_x, data_points, alpha=0.7, color='red', s=10, label=None)

        # Add labels and title
        plt.xlabel('Sequence ID (qseqid)')
        plt.ylabel('Alignment Length (bp)')
        plt.title('Box plot of Alignment Length grouped by qseqid')

        # Rotate x-axis labels if necessary (optional, depending on how many categories you have)
        plt.xticks(rotation=90)
        plt.tight_layout()
        
        output = f"{self.assay}_alignment_boxplot.png" if self.output is None else f"{self.output}_alignment_boxplot.png"
        output = os.path.join(self.dir, output)

        # Save the figure
        plt.savefig(output, dpi=120)
        
    @staticmethod
    def three_prime_mismatch(row):
        '''The Effect of Single Mismatches on Primer Extension 2018'''
        qseq_3prime_end = row['qseq'][-5:]  # Last 5 bps for forward direction
        sseq_3prime_end = row['sseq'][-5:]
        return any(q != s for q, s in zip(qseq_3prime_end, sseq_3prime_end))
    
    @staticmethod
    def pairwise_mismatch(row):
        alignment = ""
        for i, j in zip(row['qseq'], row['sseq']):
            if i == j: 
                alignment += "."
            else:
                alignment += j
        return alignment

    @staticmethod
    def proportion_test_alignment(row, snps):
        
        if row['length'] == row['qlen']:
            return pd.Series([np.inf, 1])
        
        z, p = proportions_ztest(
            count=row['length'], 
            nobs=row['qlen'], 
            value=(row['qlen'] - snps) / row['qlen'],
            #prop_var=(row['qlen'] - snps) / row['qlen'],
            alternative='smaller'
        )

        return pd.Series([round(z, 3), format(p, ".2e")])
    
    def add_primer_metadata(self):
        logging.info(f"Adding primer metadata to filtered BLAST results.")
        self.filtered_blast['3prime_mismatch'] = self.filtered_blast.apply(self.three_prime_mismatch, axis=1)
        self.filtered_blast['primer_name'] = self.filtered_blast['qseqid'].str.replace(r'_oligo\d+', '', regex=True)
        self.filtered_blast['primer_orientation'] = self.filtered_blast['primer_name'].map(self.primer_name_orientation_dict)
        self.filtered_blast['primer_assay'] = self.filtered_blast['primer_name'].map(self.primer_name_assay_dict)
        self.filtered_blast['alignment_mismatch'] = self.filtered_blast.apply(self.pairwise_mismatch, axis=1)
        self.filtered_blast['assay_type'] = 'Inclusive' if self.assay == 'i' else 'Exclusive'
        
    def check_duplicates(self):
        return self.filtered_blast.duplicated(subset=['qseqid', 'sseqid']).any()

    def save_filtered_blast_data(self):
        '''Save the filtered data to a file'''

        output = f"{self.assay}_blast_filtered_results.tsv" if self.output is None else f"{self.output}_blast_filtered_results.tsv"

        self.filtered_blast = self.filtered_blast.sort_values(by=['sseqid', 'qseqid'], ascending=[True, True])

        try:
            self.filtered_blast.to_csv(os.path.join(self.dir, output), sep='\t', index=False)
            logging.info(f"Filtered data saved to {output}.")

        except Exception as e:
            logging.error(f"Failed to save filtered data: {str(e)}")
            raise RuntimeError(f"Failed to save filtered data: {str(e)}")

    def blast_filter_workflow(self):
        logging.info(f"Filtering results for {'Inclusive' if self.assay == 'i' else 'Exclusive'} assay.")
        
        filters = {
        'seqid': self.filter_by_seqid,
        'coverage': self.filter_by_coverage_breadth,
        'uncovered': self.filter_by_coverage_uncovered,
        'mismatch': self.filter_by_mismatch,
        'snps': self.filter_by_SNPs,
        }

        logging.info(f"Number of records in raw BLAST file ({os.path.basename(self.blast_file)}) is: {self.raw_blast_data.shape[0]}")
        # Apply filters for each item
        for param, filter_func in filters.items():
            if getattr(self, param) is not None:
                filter_func()

        self.filter_highest_bitscore_per_qseqid()
        self.add_primer_metadata()
        
        if self.check_duplicates():
            logging.warning(f"Duplicate rows after filtering were found.")
        else:
            logging.info(f"No duplicate rows were found after filtering.")
        self.save_filtered_blast_data()

class BlastStatistics:
    def __init__(
            self,
            inblast: pd.DataFrame,
            bed_df: pd.DataFrame,
            oligo_dict: dict,
            reference_genomes_list: int,
            output: str,
            outdir: str,
            assay: str
        ):
        
        self.filtered_blast = inblast
        self.bed_df = bed_df
        self.oligo_dict = oligo_dict
        self.reference_genomes_list = reference_genomes_list
        self.output = output
        self.outdir = outdir
        self.assay = assay
        self.oligo_stats_df = pd.DataFrame()
        self.overall_assay_results = pd.DataFrame()
        self.primer_set_results = pd.DataFrame()

    
    def calculate_individual_oligo_stats(self):
        logging.info(f"Calculating statistics for each individual oligo.")
        try:
            total_reference_sseqid = len(self.reference_genomes_list)
            self.oligo_stats_df = self.filtered_blast.groupby(['qseqid'])['sseqid'].nunique().reset_index()
            self.oligo_stats_df.rename(columns={'qseqid': 'oligo_id', "sseqid": "pass_count"}, inplace=True)
            self.oligo_stats_df['pass_percent'] = round((self.oligo_stats_df['pass_count'] / total_reference_sseqid) * 100, 3)
            self.oligo_stats_df['total_reference'] = total_reference_sseqid
            
        except Exception as e:
            logging.error(f"An error occurred while calculating oligo statistics: {e}")
            raise RuntimeError(f"An error occurred while calculating oligo statistics: {e}")

    def calculate_overall_assay_stats(self):
        '''calculate the overall rate per assay'''
        logging.info(f"Calculating overall assay statistics.")
        
        total_reference_sseqid = len(self.reference_genomes_list)

        all_assay_data = []
        primer_assays = self.bed_df['primer_assay'].unique().tolist()
        
        for primer_assay in primer_assays:
            
            primer_assay_data = None

            # if df is empty then there are no hits for an assay
            df_primerset = self.filtered_blast[self.filtered_blast['primer_assay'] == primer_assay]

            if df_primerset.empty:
                primer_assay_data = [primer_assay, 0, 100, 0, total_reference_sseqid, total_reference_sseqid]
                all_assay_data.append(primer_assay_data)
                continue

            # check to make sure that all expected primers are found using their orientation symbol based on sseqid
            # this cannot handle duplicate hits between the qseqid found at different regions in a sseqid
            # duplicates will be removed and not considered a hit
            # however this should be improved on later with addition metrics for determining correct hit region wrt primer assay set
            primer_orientations = df_primerset.groupby(['sseqid'])['primer_orientation'].agg(list)
            orientation_elements = ['+', '.', '-']
            full_primer_assays = primer_orientations[primer_orientations.apply(lambda x: sorted(x) == sorted(list(orientation_elements)))].reset_index()
            
            if full_primer_assays.empty:
                primer_assay_data = [primer_assay, 0, 100, 0, total_reference_sseqid, total_reference_sseqid]
                all_assay_data.append(primer_assay_data)
                continue
            
            # remove those rows without all primer sequences on a sseqid
            sseqid_list = full_primer_assays['sseqid'].unique().tolist()
            df_primerset = df_primerset[df_primerset['sseqid'].isin(sseqid_list)]

            # if sstrand is minus, need to reverse the order of the sstart and ssend values for correct orientation
            # sstrand == 'minus' means that the primer was detected on the reverse compliment sseqid
            mask = df_primerset['sstrand'] == 'minus'
            df_primerset.loc[mask, ['sstart', 'send']] = df_primerset.loc[mask, ['send', 'sstart']].values

            # check primer orientation is in correct order
            ptable_sstart = pd.pivot_table(index='sseqid', columns='primer_orientation', values='sstart', fill_value=0, data=df_primerset)
            ptable_sstart = ptable_sstart.rename(columns={'+': 'sstart_+', '.': 'sstart_.', '-': 'sstart_-'})
            
            ptable_ssend = pd.pivot_table(index='sseqid', columns='primer_orientation', values='send', fill_value=0, data=df_primerset)
            ptable_ssend = ptable_ssend.rename(columns={'+': 'send_+', '.': 'send_.', '-': 'send_-'})
            
            ptable_orientation = pd.merge(ptable_sstart, ptable_ssend, on='sseqid', how='outer')
            ptable_orientation = ptable_orientation[['sstart_+', 'send_+', 'sstart_.', 'send_.', 'sstart_-', 'send_-']]
            
            def check_orientation(row):
                return row['sstart_+'] < row['send_+'] < row['sstart_.'] < row['send_.'] < row['sstart_-'] < row['send_-']
            
            ptable_orientation['in_order'] = ptable_orientation.apply(check_orientation, axis=1)
            ptable_orientation = ptable_orientation[ptable_orientation['in_order']]

            if ptable_orientation.empty:
                primer_assay_data = [primer_assay, 0, 100, 0, total_reference_sseqid, total_reference_sseqid]
                all_assay_data.append(primer_assay_data)
                continue

            # calculate statistics below after filtering above
            sseqid_list = ptable_orientation.index

            primer_assay_pass_count = len(sseqid_list)
            primer_assay_fail_count = total_reference_sseqid - primer_assay_pass_count
            primer_filtered_pass_percent = round(primer_assay_pass_count / total_reference_sseqid * 100, 3)
            primer_filtered_fail_percent = round((100 - primer_filtered_pass_percent), 3)
            
            primer_assay_data = [
                primer_assay,
                primer_filtered_pass_percent,
                primer_filtered_fail_percent,
                primer_assay_pass_count,
                primer_assay_fail_count,
                total_reference_sseqid
            ]

            all_assay_data.append(primer_assay_data)

            try:
                self.primer_clustermap(df=df_primerset, primer_assay=primer_assay)
                logging.info(f"Cluster map created.")
            except ValueError:
                logging.error(f"Unable to create clustermap.")
                pass

        columns = [
            "primer_assay",
            "assay_pass_percent",
            "assay_fail_percent",
            "assay_pass_count",
            "assay_fail_count",
            "total_references"
        ]

            
        
        self.overall_assay_results = pd.DataFrame(all_assay_data, columns=columns).sort_values(by=['primer_assay'])
        
    def calculate_oligo_set_stats(self):
        logging.info(f"Calculating statistics for oligo sets.")
        
        data = {
            'primer_assay':[],
            'oligo_set': [], 
            'pass_percent': [],
            'fail_percent': [],
            'pass_count': [],
            'fail_count': [],
            'filtered_reference_count': [],
            'forward': [], 
            'reverse': [], 
            'probe': []
        }
        
        total_reference_sseqid = len(self.reference_genomes_list)
        all_assay_data = []
        primer_assays = self.bed_df['primer_assay'].unique().tolist()
    
        for primer_assay in primer_assays:
            
            primer_assay_data = None

            # if df is empty then there are no hits for an assay
            df_primerset = self.filtered_blast[self.filtered_blast['primer_assay'] == primer_assay]

            if df_primerset.empty:
                primer_assay_data = [primer_assay, pd.NA, 0, 100, 0, total_reference_sseqid, total_reference_sseqid, pd.NA, pd.NA, pd.NA]
                all_assay_data.append(primer_assay_data)
                continue
            
            forward_list = self.oligo_dict[primer_assay]['Forward']
            probe_list = self.oligo_dict[primer_assay]['Probe']
            reverse_list = self.oligo_dict[primer_assay]['Reverse']

            oligo_primer_sets = list(product(forward_list, probe_list, reverse_list))
            
            for idx, oligo_set in enumerate(oligo_primer_sets, start=1):

                oligo_set_idx = f"{primer_assay}_primerset_{idx}"
                # check to make sure that all expected primers are found using their orientation symbol based on sseqid
                # this cannot handle duplicate hits between the qseqid found at different regions in a sseqid
                # duplicates will be removed and not considered a hit
                # however this should be improved on later with addition metrics for determining correct hit region wrt primer assay set
                
                # keep only oligos in a set per primer assay
                df_primerset = df_primerset[df_primerset['qseqid'].isin(oligo_set)]
                primer_orientations = df_primerset.groupby(['sseqid'])['primer_orientation'].agg(list)
                orientation_elements = ['+', '.', '-']
                full_primer_assays = primer_orientations[primer_orientations.apply(lambda x: sorted(x) == sorted(list(orientation_elements)))].reset_index()
                
                if full_primer_assays.empty:
                    primer_assay_data = [primer_assay, oligo_set_idx, 0, 100, 0, total_reference_sseqid, total_reference_sseqid, oligo_set[0], oligo_set[1], oligo_set[2]]
                    all_assay_data.append(primer_assay_data)
                    continue
                
                # remove those rows without all primer sequences on a sseqid
                sseqid_list = full_primer_assays['sseqid'].unique().tolist()
                df_primerset = df_primerset[df_primerset['sseqid'].isin(sseqid_list)]

                # if sstrand is minus, need to reverse the order of the sstart and ssend values for correct orientation
                # sstrand == 'minus' means that the primer was detected on the reverse compliment sseqid
                mask = df_primerset['sstrand'] == 'minus'
                df_primerset.loc[mask, ['sstart', 'send']] = df_primerset.loc[mask, ['send', 'sstart']].values

                # check primer orientation is in correct order
                ptable_sstart = pd.pivot_table(index='sseqid', columns='primer_orientation', values='sstart', fill_value=0, data=df_primerset)
                ptable_sstart = ptable_sstart.rename(columns={'+': 'sstart_+', '.': 'sstart_.', '-': 'sstart_-'})
                
                ptable_ssend = pd.pivot_table(index='sseqid', columns='primer_orientation', values='send', fill_value=0, data=df_primerset)
                ptable_ssend = ptable_ssend.rename(columns={'+': 'send_+', '.': 'send_.', '-': 'send_-'})
                
                ptable_orientation = pd.merge(ptable_sstart, ptable_ssend, on='sseqid', how='outer')
                ptable_orientation = ptable_orientation[['sstart_+', 'send_+', 'sstart_.', 'send_.', 'sstart_-', 'send_-']]
                
                def check_orientation(row):
                    return row['sstart_+'] < row['send_+'] < row['sstart_.'] < row['send_.'] < row['sstart_-'] < row['send_-']
                
                ptable_orientation['in_order'] = ptable_orientation.apply(check_orientation, axis=1)
                ptable_orientation = ptable_orientation[ptable_orientation['in_order']]

                if ptable_orientation.empty:
                    primer_assay_data = [primer_assay, oligo_set_idx, 0, 100, 0, total_reference_sseqid, total_reference_sseqid, oligo_set[0], oligo_set[1], oligo_set[2]]
                    all_assay_data.append(primer_assay_data)
                    continue

                # calculate statistics below after filtering above
                sseqid_list = ptable_orientation.index

                primer_assay_pass_count = len(sseqid_list)
                primer_assay_fail_count = total_reference_sseqid - primer_assay_pass_count
                primer_filtered_pass_percent = round(primer_assay_pass_count / total_reference_sseqid * 100, 3)
                primer_filtered_fail_percent = round((100 - primer_filtered_pass_percent), 3)
                
                primer_assay_data = [
                    primer_assay,
                    oligo_set_idx,
                    primer_filtered_pass_percent,
                    primer_filtered_fail_percent,
                    primer_assay_pass_count,
                    primer_assay_fail_count,
                    total_reference_sseqid,
                    oligo_set[0],
                    oligo_set[1],
                    oligo_set[2]

                ]

                all_assay_data.append(primer_assay_data)
            
        columns = [
            "primer_assay",
            "primer_set_idx",
            "assay_pass_percent",
            "assay_fail_percent",
            "assay_pass_count",
            "assay_fail_count",
            "total_references",
            "forwared",
            "probe",
            "reverse"
        ]

        self.primer_set_results = pd.DataFrame(all_assay_data, columns=columns)


    def filter_min_bitscore(self) -> pd.DataFrame:
        '''
        only get the minimal matching score for the reduce primer matches between qseqid and sseqid
        this will not remove duplicates if there are identical bitscores.
        can use this to identify primers aligning to mutiple locations in a genome
        '''
        idxmin_ = self.filtered_blast.groupby(['qseqid', 'sseqid'])['bitscore'].idxmin()
        self.filtered_blast = self.filtered_blast.loc[idxmin_].reset_index(drop=True)
        logging.info(f"Filtered for lowest bitscore for 'qseqid' and 'sseqid'. Current records now {len(self.filtered_blast)}.")

    def filter_max_bitscore(self) -> pd.DataFrame:
        '''
        only get the maximum matching score for the reduce primer matches bettween qseqid and sseqid
        this will not remove duplicates if there are identical bitscores.
        can use this to identify primers aligning to mutiple locations in a genome
        '''
        idxmax_ = self.filtered_blast.groupby(['qseqid', 'sseqid'])['bitscore'].idxmax()
        self.filtered_blast = self.filtered_blast.loc[idxmax_].reset_index(drop=True)
        logging.info(f"Filtered for highest bitscore for 'qseqid' and 'sseqid'. Current records now {len(self.filtered_blast)}.")

    def primer_clustermap(self, df:pd.DataFrame, primer_assay:str):

        sys.setrecursionlimit(10000)
        logging.getLogger('matplotlib').setLevel(logging.CRITICAL)


        ptable = pd.pivot_table(data=df, index="sseqid", columns='primer_name', values='pident', fill_value=0)
        mask = ptable == 0
        #cmap = plt.cm.mako
        cmap = sns.color_palette("rocket_r", as_cmap=True)
        cmap.set_bad('grey')
    
        g = sns.clustermap(
            ptable,
            row_cluster=True,
            col_cluster=False,
            mask=mask,
            cmap=cmap, 
            cbar_kws={
                "shrink": 1
            }
        )

        g.ax_heatmap.set_xlabel(f"Query Sequences (qseqid, n={len(ptable.columns)})", fontsize=12)  # Change the x-axis label
        g.ax_heatmap.set_ylabel(f"Subject Sequences (sseqid, n={len(ptable.index)})", fontsize=12)  # Change the y-axis label
        
        # remove y labels if there are > 1000 genomes
        if len(ptable.index) > 1000:
            g.ax_heatmap.set_yticklabels([])
        
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        plt.tight_layout()
        output = os.path.join(self.outdir, f"{primer_assay}_clustermap.png")
        plt.savefig(output, dpi=300)

    def save_primer_set_stats_data(self):
        extention = "primer_set_statistics.tsv"
        output = f"{self.assay}_{extention}" if self.output is None else f"{self.output}_{extention}"
        pathway = os.path.join(self.outdir, output)
        self.primer_set_results.to_csv(pathway, sep='\t', index=False)
        logging.info(f"Saving Primer Set statistics to {output}")

    def save_primer_overall_stats_data(self):
        extention = "primer_overall_statistics.tsv"
        output = f"{self.assay}_{extention}" if self.output is None else f"{self.output}_{extention}"
        pathway = os.path.join(self.outdir, output)
        self.overall_assay_results.to_csv(pathway, sep='\t', index=False)
        logging.info(f"Saving Overall statistics to {output}")

    def save_primer_individual_stats_data(self):
        extention = "primer_oligo_statistics.tsv"
        output = f"{self.assay}_{extention}" if self.output is None else f"{self.output}_{extention}"
        pathway = os.path.join(self.outdir, output)
        self.oligo_stats_df.to_csv(pathway, sep='\t', index=False)
        logging.info(f"Saving Oligo statistics to {output}")

def parse_args():
    """
    Parse command-line arguments.

    :return: Namespace with parsed arguments.
    """
    prog = "IE_assay.py"
    description = \
    f'''
    This script is designed for conducting in silico Inclusivity/Exclusivity assays.  
    It handles nucleotide primer sequence disambiguation, performs local alignments using BLASTn, filters results based on specified thresholds, and calculates performance statistics of the primers.
    Note: Currently there are not paramters to change the BLAST paramters.
    ''' 
    epilog = '''
    Minimum Usage:
        IE_assay.py -1 <bed-like_primer_file.tsv> -2 <database.fasta> -a {i,e}
    '''

    parser = argparse.ArgumentParser(
        description=description,
        prog=prog,
        epilog=epilog)

    # Required arguments
    required_group = parser.add_argument_group('Required Files')
    required_group.add_argument(
        '-1', '--primers',
        type=str,
        required=True,
        help='Pathway to the primer (bed-like) sequence file and metadata.',
        dest='primer_seqs'
    )
    required_group.add_argument(
        '-2', '--database',
        type=str,
        required=True,
        help='Pathway to genomes fasta file.',
        dest='database'
    )
    required_group.add_argument(
        '-a', '--assay',
        type=str,
        required=True,
        # default='i',
        choices=['i', 'e'],
        help='Inclusive (i) or Exclusive (e) assay [-a i OR -a e]. Default parameters below for each assay.',
        dest='assay'
    )

    assay_group = parser.add_argument_group('Assay Type')
    
    assay_group.add_argument(
        '-i', '--seqid',
        help="Minimum percent sequnece ideneity between 0-100 [default: i=None, e=80]",
        default=None,
        type=float,
        dest="seqid",
        required=False,
    )

    assay_group.add_argument(
        '-c', '--coverage',
        help="Minimum sequencing breadth of coverage (length/qlen). Value between 0-1  [default: i=None, e=0.8]",
        default=None,
        type=float,
        dest="coverage",
        required=False,
    )

    assay_group.add_argument(
        '-u', '--uncovered',
        help="Maximum number of uncovered bases [default: i=2, e=None]",
        default=None,
        required=False,
        type=int,
        dest="uncovered"
    )

    assay_group.add_argument(
        '-m', '--mismatch',
        help="Maximum number of mismatch [default: i=2, e=None]",
        default=None,
        type=int,
        dest="mismatch",
        required=False,
    )

    assay_group.add_argument(
        '-s', '--snps',
        help="Maximum number of snps (mismatch + gaps) [default: i=2, e=None]",
        type=int,
        dest="snps",
        required=False,
        default=None
    )

    # Additional arguments
    optional_group = parser.add_argument_group('Additional arguments')

    optional_group.add_argument(
        '-d', '--dir',
        type=str,
        required=False,
        default="xSILICO_RESULTS",
        help='Name of directory for storing data (default: xSILICO_RESULTS): NOTE: Will Overwrite any preexisting Directory!!! ***User Beware***',
        dest='dir'
    )

    optional_group.add_argument(
        '-p', '--pathway',
        type=str,
        required=False,
        default=os.getcwd(),
        help='Pathway to save storage directory and results (default: current working directory)',
        dest='pathway'
    )

    optional_group.add_argument(
        '-o', '--output',
        type=str,
        required=False,
        default=None,
        help='Prefix given to all BLAST results.',
        dest='output'
    )

    optional_group.add_argument(
        '-t', '--threads',
        type=str,
        default='1',
        help='Number of threads for blastn [default: 1]',
        dest='threads'
    )

    optional_group.add_argument(
        '-l', '--logname',
        type=str,
        default=None,
        help='Prefix name for logging file (default: IE_assay.log)',
        dest='logname'
    )

    optional_group.add_argument(
        '-C', '--check',
        dest='check_fasta',
        action='store_true',
        help='Check ALL records in fasta files are formated correctly (NOTE: may fail due to memory issue with larger fasta files)'
    )

    optional_group.add_argument(
        '-v', '--verbose',
        dest='verbose',
        action='store_true'
    )

    opts = parser.parse_args()

    return opts

def main(tempdir):

    args = parse_args()
    # required arguments
    primer_seqs = args.primer_seqs
    database = args.database
    assay = args.assay
    # assay arguments
    seqid = args.seqid
    coverage = args.coverage
    uncovered = args.uncovered
    mismatch = args.mismatch
    snps = args.snps
    # optional arguments
    threads = args.threads
    parent_dir = args.dir 
    pathway = args.pathway
    logname = args.logname
    output = args.output
    check_fasta = args.check_fasta
    verbose = args.verbose

    assay_type = {'i': "Inclusivity", 'e': "Exclusivity"}
    logname = setup_logging_info(logname, verbose=verbose)
    logging.info(f"Starting {assay_type[assay]} Run")
    logging_args(args)

    if check_fasta:
        logging.info("Checking Input Files...")
        check_file_existance(primer_seqs)
        check_file_existance(database)
        check_fasta_existance(database)
        logging.info(f"Input files passed check points.")
    else:
        logging.info("Bypassing fasta file check.")
    

    logging.info(f"Disambiguating Degenerate Sequences.")
    dsp = DegenerateSequenceProcessor(
        bedfile=primer_seqs, 
        dir=temp_dir
    )
    dsp.bed2df()
    dsp.bed2seqs()
    dsp.process_degenerate_sequences()
    dsp.write_degenerate_sequences()

    logging.info(f"Running BLAST Utilities.")  
    blast_util = BlastUtilities(
        query=dsp.output_fasta,
        reference=database,
        dir=temp_dir,
        assay=assay,
        output=output,
        threads=threads
    )

    blast_util.is_blastn_available()
    blast_util.make_blast_db()
    blast_util.run_blastn()


    logging.info(f"Filtering raw BLAST Data.")
    blast_filter = BlastFilter(
        blast_file=blast_util.blast_file, # call raw blast file from blastuntilities 
        output=output,
        dir=temp_dir,
        primer_name_assay_dict = dsp.primer_name_assay_dict(),
        primer_name_orientation_dict = dsp.primer_name_orientation_dict(),
        assay=assay,
        seqid=seqid,
        coverage=coverage,
        mismatch=mismatch,
        uncovered=uncovered,
        snps=snps
        )
    
    blast_filter.load_assay_parameters()
    blast_filter.load_data()
    blast_filter.alignment_distribution()
    # self.alignment_boxplot()          # workshop later
    blast_filter.blast_filter_workflow()

    logging.info(f"Calculating statistics from BLAST filtered data.")

    blast_stats = BlastStatistics(
        inblast=blast_filter.filtered_blast,
        bed_df=dsp.bed_df, #bed-like file for sotring information
        oligo_dict=dsp.oligo_dict,
        reference_genomes_list=blast_util.extract_reference_deflines(),
        output=output,
        outdir=temp_dir,
        assay=assay
    )

    blast_stats.calculate_individual_oligo_stats()
    # create new column with primer orientation and primer name for overall statistics
    blast_stats.filter_max_bitscore() 

    blast_stats.calculate_overall_assay_stats()
    blast_stats.calculate_oligo_set_stats()

    blast_stats.save_primer_overall_stats_data()
    blast_stats.save_primer_set_stats_data()
    blast_stats.save_primer_individual_stats_data()
    #blast_stats.primer_clustermap()

    logging.info(f"Finished.")

    target_dir = os.path.join(pathway, parent_dir)
    if os.path.isdir(target_dir):
        # Remove the existing target directory if it exists
        shutil.rmtree(target_dir)
    # Move the temporary directory to the target location

    shutil.move(temp_dir, target_dir)
    shutil.move(logname, target_dir)
    

if __name__ == "__main__":

    random_chars = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    with tempfile.TemporaryDirectory(suffix='_xsilicotemp', prefix=random_chars, dir='.') as temp_dir:
        main(temp_dir)