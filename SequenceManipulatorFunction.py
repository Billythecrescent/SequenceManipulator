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
from Bio.SeqUtils import GC
from Bio.Alphabet import *
import pandas as pd
import numpy as np

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

def fastaSeqsAttribute(path):
    '''Find the IDs, lengths and GC contents of all the sequences in a fasta file.
    path: string
        The path of the fasta file

    Return:
    ------
    combined: DataFrame
        The attributes of all the sequences of the fasta file.
    '''
    IDs = []
    lengths = []
    descriptions = []
    GCs = []
    with open(os.path.join(path), "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            IDs.append(record.id)
            lengths.append(len(record.seq))
            GCs.append(round(GC(record.seq), 2))
            descriptions.append(record.description)
    
    combined = np.array([IDs, lengths, GCs, descriptions])
    combined = combined.T
    df = pd.DataFrame(combined, columns=['sequenceID', 'length', 'GC', 'description'])
    return df

def fastaSpecifyAlphabet(path, alphabet):
    '''Read and specify the alphabet of all the sequences in a fasta file.
    path: string
        The path of the fasta file
    alphabet: Bio.Alphabet
        The alphabet to be specified to the sequences.
    
    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The complemented sequences of all DNA sequences in a fasta file.
    '''

def fastaComplement(path, out=None):
    '''Find the complements of DNA sequences in a fasta file.
    path: string
        The path of the fasta file
    out: string
        output filename, containing '.fasta' or not is both acceptable
        output format: FASTA

    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The complemented sequences of all DNA sequences in a fasta file.
    '''
    seqs = readFASTA(path)
    for key in seqs:
        seqs[key] = seqs[key].complement()
    if out != None:
        outFile = out if out.split('.')[-1] == 'fasta' else out + ".fasta"
        SeqIO.write(seqs, outFile, "fasta")
    return seqs

def fastaReverse(path, out=None):
    '''Find the reverse complements of DNA sequences in a fasta file.
    path: string
        The path of the fasta file
    out: string
        output filename, containing '.fasta' or not is both acceptable
        output format: FASTA

    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The reverse-complemented sequences of all DNA sequences in a fasta file.
    '''
    seqs = readFASTA(path)
    for key in seqs:
        seqs[key] = seqs[key][::-1]
    
    if out != None:
        outFile = out if out.split('.')[-1] == 'fasta' else out + ".fasta"
        SeqIO.write(seqs, outFile, "fasta")
    return seqs

def fastaReverseComplement(path, out=None):
    '''Find the reverse of DNA sequences in a fasta file.
    path: string
        The path of the fasta file
    out: string
        output filename, containing '.fasta' or not is both acceptable
        output format: FASTA

    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The reversed sequences of all DNA sequences in a fasta file.
    '''
    seqs = readFASTA(path)
    for key in seqs:
        seqs[key] = seqs[key].reverse_complement()
    
    if out != None:
        outFile = out if out.split('.')[-1] == 'fasta' else out + ".fasta"
        SeqIO.write(seqs, outFile, "fasta")
    return seqs

print(fastaSeqsAttribute(os.path.join(seq_path, "VACV_virus.fasta")))
# fastaReverseComplement(os.path.join(seq_path, "E.setchuanus_cytb.fasta"))


