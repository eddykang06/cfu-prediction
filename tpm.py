# Script to extract Fcnts data from a folder, CFU data from a file, then
# bind everything into 1 dataframe that can be used for training regression models
# 
# Input: 
#      
# Note that the input folder should be the default


# Packages
import pandas as pd
import numpy as np
import os
import glob
import sys

# Store command-line input and output
fcnts_path = sys.argv[0] # folder
cfu_path   = sys.argv[1]
outfig     = sys.argv[2]

# Function to load 


# can use df.rename(lambda x: ) to rename the samples?
# Function to load the Fcnts data
def read_fcnts(folder_path):

    # Extract all file names within folder
    files = 


# Function to load the CFU data
def read_cfu(file_path):

# Function to convert Fcnt files to dataframes of TPMs
def tpm_covert(fcnt_files):


# Function to bind TPMs
def bind_data():

# Call main
def main(fcnts_path, cfu_path, outfig):
