'''
File: SeqMain.py
Author: Mucun Hou
Date: May 21, 2020
Description: This script is a main script used to build the interface between 
    the SequenceManipulatorFunction and the user
'''

from Bio.Seq import Seq
import os.path
import SequenceManipulatorFunction as SMF

##--- File Paths ---##
module_path = os.path.abspath(os.path.dirname(__file__)) #code\MHC-peptide_prediction
seq_path = os.path.join(module_path, "sequences")

def main(ifFile, inputStr, outputFile, index=-1, inputFormat='fasta', outputFormat='fasta'):

    #Intro
    print()
    print('***************************************')
    print('*   Welcome to Sequence Manipulator   *')
    print('***************************************')

    SingleFlag = 0

    if not ifFile:
        print("*   Your input is a sequence   *\nSeq: %s" %(inputStr + "..." if len(inputStr) > 100 else inputStr))
        SingleFlag = 1
    else:
        print("*   Your input is a file   *\nPath: %s" %inputStr)
        sequences = SMF.readFASTA(inputStr)
        if index < -1:
            print("ERROR: index not valid! Please specify zero or positive integar.")
            return
        elif index >= 0:
            print("*   You choose sequence %d   *\nSeq: %s" %(index, (sequences[index] + "..." if len(sequences[index]) > 100 else sequences[index])))
            SingleFlag = 1
        else:
            print("*   Your choose all sequences in the file   *")    

    if SingleFlag:
        validFunction = ["SequenceLength", "fastaSeqsAttribute", "Complement", "ReverseComplement", "fastaTranscribe", "fastaReverseTrasncribe", "fastaTranslate"]

