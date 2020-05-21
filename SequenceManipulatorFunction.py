'''
File: SequenceManipulatorFunction.py
Author: Mucun Hou
Date: May 12, 2020
Description: This script is to provide some easy-to-use tools
    for sequence manipulation.
'''
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

##--- File Paths ---##
module_path = os.path.abspath(os.path.dirname(__file__)) #code\MHC-peptide_prediction
seq_path = os.path.join(module_path, "sequences")

def readFASTA(path):
    '''read the fasta file by a handler and return the list of
    Seq record objects.
    path: string
        The path of the fasta file

    Return:
    ------
    sequences: dict[Seq record]
        The dictionary of Seq record with the id as the key. 
    '''
    sequences = {}
    with open(os.path.join(path), "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = record.seq
    return sequences

def fastaNumber(path):
    '''Find how many sequences in a fasta file
    path: string
        The path of the fasta file

    Return:
    ------
    num: int
        The sequence number of the fasta file
    '''
    seqs = readFASTA(path)
    num = len(seqs)
    return num

def fastaLength(path):
    '''File the whole length of all the sequences in a fasta file.
    path: string
        The path of the fasta file

    Return:
    ------
    length: int
        The overall length of all sequences of the fasta file.
    '''
    seqs = readFASTA(path)
    length = 0
    for key in seqs:
        length = length + len(seqs[key])
    return length

