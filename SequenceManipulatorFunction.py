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
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_nucleotide
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

def fastaSpecifyAlphabet(path, alphabet=None):
    '''Read and specify the alphabet of all the sequences in a fasta file.
    path: string
        The path of the fasta file
    alphabet: Bio.Alphabet.IUPAC
        The IUPAC alphabet to be specified to the sequences.

    RUN fastaSpecifyAlphabet(path,alphabet=None) for available alphabet help.

    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The complemented sequences of all DNA sequences in a fasta file.
    '''
    if alphabet == None:
        print("Available alphabets: IUPAC.ambiguous_dna, IUPAC.unambiguous_dna, IUPAC.ambiguous_rna, IUPAC.unambiguous_rna, IUPAC.protein, IUPAC.extended_protein, generic_nucleotide, SingleLetterAlphabet, NucleotideAlphabet, ProteinAlphabet, DNAAlphabet, RNAAlphabet.")
        return
    sequences = {}
    with open(os.path.join(path), "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq
            try:
                seq.alphabet = alphabet
            except AttributeError:
                
                return
            else:
                sequences[record.id] = seq
    print(sequences)
    return sequences

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

def fastaTranscribe(path, transcribe_options=None, out=None):
    '''Transcribe all the sequences in a fasta file. Those are seemed as coding sequences
    if no transcribe_options are given. Users can specify whether the sequences are coding
    sequences or template sequences by transcribe_options. 'c' or 'C' means coding sequence,
    and 't' or 'T' means template sequence. For example, ['C', 'C', 'T'] means that there 
    are three sequences and first two are coding sequences and the last is a template sequence.

    path: string
        The path of the fasta file
    transcribe_options: list or tuple
        Specify the sequences in the file as coding sequence ('C' or 'c') or template sequence
        ('T' or 't'). ['C', 'C', 'T'] for example.
    out: string
        output filename, containing '.fasta' or not is both acceptable
        output format: FASTA

    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The transcribed sequences of all RNA sequences in a fasta file.
    '''
    seqs = readFASTA(path)
    RNAs = {}
    if transcribe_options == None:
        for key in seqs:
            rna = seqs[key].transcribe()
            RNAs[key] = rna
    elif type(transcribe_options) != list and type(transcribe_options) != tuple:
        print("ERROR: Please specify the transcribe_options by list or tuple.")
        return
    elif len(transcribe_options) != len(seqs):
        print("ERROR: The length of transcribe_options %d is not consistent with the number of sequences %d" %(len(transcribe_options), len(seqs)))
        return
    else:
        with open(os.path.join(path), "rU") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                key = record.id
                seq = record.seq
                if transcribe_options[i] == 'T' or transcribe_options[i] == 't':
                    rna = seq.reverse_complement().transcribe()
                elif transcribe_options[i] == 'C' or transcribe_options[i] == 'c':
                    rna = seq.transcribe()
                else:
                    print("ERROR: Illegal character %c in transcribe_options list" %transcribe_options[i])
                    return
                RNAs[key] = rna
    #ouput FASTA file
    if out != None:
        outFile = out if out.split('.')[-1] == 'fasta' else out + ".fasta"
        SeqIO.write(RNAs, outFile, "fasta")

    return RNAs

def fastaReverseTrasncribe(path, transcribe_options=None, out=None):
    '''Reverse-transcribe all the sequences in a fasta file. Sequences can be reverse-transcribed 
    to coding sequences or template sequences by specifying the transcribe-options. If no given, 
    this function will convert the RNAs to coding sequences. In transcribe_options. 'c' or 'C' 
    means coding sequence, and 't' or 'T' means template sequence. For example, ['C', 'C', 'T'] 
    means that there are three sequences and first two are to coding sequences and the last is to 
    a template sequence.

    path: string
        The path of the fasta file
    transcribe_options: list or tuple
        Specify the RNA sequences the RNA sequences in the file are reverse-transcribed to.
        Whether to coding sequence ('C' or 'c') or to template sequence ('T' or 't'). 
        ['C', 'C', 'T'] for example.
    out: string
        output filename, containing '.fasta' or not is both acceptable
        output format: FASTA

    Return:
    ------
    seqs: dictionary{seqID: Seq}
        The transcribed sequences of all DNA sequences in a fasta file.
    '''
    seqs = readFASTA(path)
    DNAs = {}
    if transcribe_options == None:
        for key in seqs:
            dna = seqs[key].back_transcribe()
            DNAs[key] = dna
    elif type(transcribe_options) != list and type(transcribe_options) != tuple:
        print("ERROR: Please specify the transcribe_options by list or tuple.")
        return
    elif len(transcribe_options) != len(seqs):
        print("ERROR: The length of transcribe_options %d is not consistent with the number of sequences %d" %(len(transcribe_options), len(seqs)))
        return
    else:
        with open(os.path.join(path), "rU") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                key = record.id
                seq = record.seq
                if transcribe_options[i] == 'T' or transcribe_options[i] == 't':
                    dna = seq.back_transcribe().reverse_complement()
                elif transcribe_options[i] == 'C' or transcribe_options[i] == 'c':
                    dna = seq.back_transcribe()
                else:
                    print("ERROR: Illegal character %c in transcribe_options list" %transcribe_options[i])
                    return
                DNAs[key] = dna
    #ouput FASTA file
    if out != None:
        outFile = out if out.split('.')[-1] == 'fasta' else out + ".fasta"
        SeqIO.write(DNAs, outFile, "fasta")

    print(DNAs)
    return DNAs

# print(fastaSeqsAttribute(os.path.join(seq_path, "VACV_virus.fasta")))
# print(fastaNumber(os.path.join(seq_path, "E.setchuanus_cytb.fasta")))
# fastaTranscribe(os.path.join(seq_path, "E.setchuanus_cytb.fasta"), transcribe_options=['T']*27)


