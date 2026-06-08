"""All data extraction functions for models using DEG data to predict synergy"""
import pandas as pd
import numpy as np
import os
import functools
from functools import reduce

def clean_csv_name(name):
    """
    Function to clean a csv name into a readable sample ID

    Args:
        name : ex: "T4-wt1AMX1hr-T4-wtNDC0hr-DGE-results"
    
    Returns:
        cleaned : ex: "1AMX1hr"
    """
    suffixes = [
        ".csv",
        "-T4-wtNDC0hr-DGE-results",
        "-T4-wtNDC1hr-DGE-results",
        "-T4-wtNDC2hr-DGE-results",
        "-T4-wtNDC4hr-DGE-results",
        "T4-wt"
    ]
    
    cleaned = name

    for suffix in suffixes:
        cleaned = cleaned.replace(suffix, "")

    return cleaned

def read_l2fc_as_df(data_dir, time_matched):
    """
    Function to extract all l2fc values from a directory of DGE results for various conditions.
    The output will be a list of dataframes, where each DGE result is stored as a dataframe with gene IDs on index, and 1 column of l2fc values

    Args:
        data_dir     : path to folder containing DGE results
        time_matched : boolean indicating whether DGE results should be from time zero or time-matched NDC

    Returns:
        l2fc_df_list : list of dataframes storing l2fc results
        names        : sample IDs corresponding to each dataframe
    """
    all_files = os.listdir(data_dir)
    files = []

    # Time-matched NDC comparisons
    if time_matched:
        for f in all_files:
            if f.endswith(".csv") and "NDC0hr" in f:
                files.append(f)\
    
    # Time zero comparisons
    else:
        for f in all_files:
            if f.endswith(".csv") and "NDC0hr" not in f:
                files.append(f)

    if len(files) == 0: 
        raise FileNotFoundError(f"No matching DGE files found in {data_dir}")

    # Sort and specify path to data files
    files = sorted(files)

    # Get sample IDs
    ids = [clean_csv_name(f) for f in files]

    # Convert csvs to dataframes
    filenames = ["".join([data_dir, "/" , f]) for f in files]
    l2fc_df_list = [pd.read_table(filenames[i], 
                                  sep = ",", 
                                  header = 0, 
                                  index_col = 1)["log2FoldChange"].rename(ids[i]) # Genes
                                  for i in range(len(filenames))]

    return l2fc_df_list, ids

def bind_l2fc_data(l2fc_df_list, ids):
    """
    Function to take a list of l2fc DEG dataframes, then bind all into 1 dataframe
    Args: 
        l2fc_df_list : list of l2fc dataframes
        ids          : corresponding sample IDs

    Returns:
        all_l2fc [N,G] : dataframe with all feature values (N samples on row, G genes on column)
    """
    all_l2fc = reduce(lambda df1, df2 :
                      pd.merge(df1, df2, 
                               left_index = True, 
                               right_index = True, 
                               how = "outer"),
                        l2fc_df_list)

    # Tranpose to get genes on columns
    all_l2fc = all_l2fc.T

    return all_l2fc


def read_avg_cfus(folder_path):
    """
    Function to extract CFUs into a dataframe with conditions on rows, 1 CFU column
    Args:
        folder_path : path to folder containing CFUs (same format as /all_cfus)

    Returns:
        all_avg_cfus [N,1] : df with condition names as index, 1 column of avg CFUs across triplicates (N = # samples)
    """

    # Get files
    files = os.listdir(folder_path)
    
    # Select CSV files
    cfu_files = [csv for csv in files if ".csv" in csv]
    cfu_files = sorted(cfu_files)
    cfu_files = ["".join([folder_path, "/", csv]) for csv in files]
    
    # Load each file as a dataframe
    cfu_dfs = [pd.read_table(csv, sep = ",", header = 0) for csv in cfu_files] 

    # change all of this to calculate average CFU for each condition
    for i, df in enumerate(cfu_dfs):

        # Remove the triplicate column
        df = df.drop(columns = "Triplicates")

        # Calculate mean of three triplicates and store
        means = df.mean()
        cfu_dfs[i] = means
    
    # Concat
    all_avg_cfus = pd.concat(cfu_dfs, axis = 0) # series
    all_avg_cfus = all_avg_cfus.rename("CFU").to_frame()

    # Remove space breaker characters
    all_avg_cfus.index = [id.replace("\xa0", "") for id in all_avg_cfus.index]
    
    return all_avg_cfus


def bind_all_data(feature_df, cfu_df):
    """
    Function bind TPM and cfu dfs
    Args:
        featurem_df [N,G] : Dataframe of TPMs, N = # samples, G = # genes, labels on index
        cfu_df [N,1]      : Dataframe of CFUs, labels on index

    Returns:
        data_df : [N, G+1] : Dataframe of all TPMs and CFUs as last column
    """
    # Right join so that CFUs exist
    data_df = pd.merge(feature_df, cfu_df, left_index = True, right_index = True, how = "right")
    data_df = data_df.iloc[data_df.index.argsort()]

    return data_df


def get_l2fc_and_cfu_data(l2fc_dir, cfu_dir, time_matched):
    
    l2fc_df_list, ids = read_l2fc_as_df(l2fc_dir, time_matched)
    all_l2fc = bind_l2fc_data(l2fc_df_list, ids)
    all_avg_cfus = read_avg_cfus(cfu_dir)
    data_df = bind_all_data(all_l2fc, all_avg_cfus)
    
    return data_df


"""Functions to calcualte transcriptional interaction scores"""
def simple_interaction_score(l2fc_a, l2fc_b, l2fc_ab):
    """
    Function to a compute a simple interaction score to quantify non-additivity in gene expression 
    using single drug log2foldchange and combination log2foldchange

    Args:
        l2fc_a   : Log2foldchange in condition A
        l2fc_b   : Log2foldchange in condition B
        l2fc_ab : Log2foldchange in combination

    Returns:
        score : Simple interaction score
    """
    score = l2fc_ab - l2fc_a - l2fc_b

    return score


"""Functions to calculate drug synergy scores"""
def eob_score(v_a, v_b, v_ab):
    """
    Function to calculate the excess over bliss (EOB) synergy score for given single drug survival and a 
    combination drug survival

    Args:
        v_a  : Fraction surviving cells in condition A
        v_b  : Fraction surviving cells in condition B
        v_ab : Fraction surviving cells in condition A+B

    Returns:
        score: EOB score calculated using fraction of surviving cells
    """
    # EOB score (calculated using survival fraction)
    score = v_a * v_b - v_ab

    return score


def hsa_score(v_a, v_b, v_ab):
    """  
    Function to calculate the HSA synergy score for given single drug survival and a 
    combination drug survival

    Args:
        v_a  : Fraction surviving cells in condition A
        v_b  : Fraction surviving cells in condition B
        v_ab : Fraction surviving cells in condition A+B

    Returns:
        score: HSA score calculated using fraction of surviving cells
    """
    # HSA score (calculated using survival fraction)
    score = np.min([v_a, v_b]) - v_ab

    return score


"""Functions to generate new data df for synergy prediction"""