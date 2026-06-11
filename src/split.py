"""Functions to implement specific train-test splits"""

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

def combination_stratified_split(df, num_folds, seed = None):
    """
    Train-test split where train-test split is stratified by combination, so the 
    drug distribution is equal between train-test

    Args:
        df : Dataframe with metadata attached
    
    Returns:
        splits : List of tuples with (train_idx, test_idx)
    """
    # Store drug IDs to stratify by 
    ids = df["drug_id"].to_numpy()

    # Initialize splitter
    cv = StratifiedKFold(
        n_splits = num_folds,
        shuffle = True,
        random_state = seed
    )

    # Stratified split
    splits = list(cv.split(df, ids))

    return splits


def combination_held_out_split(df):
    """
    Train-test split where 1 drug combination is held out in each split

    Args:
        df : Dataframe with metadata attached
    
    Returns:
        splits : List of tuples with (train_idx, test_idx)
    """

    # Store drug IDs and select unique combinations
    ids = df["drug_id"].to_numpy()
    unique_combos = np.unique(ids)

    # Store splits
    splits = []

    for combo in unique_combos:
        
        # Get idx
        train_idx = np.array([i for i in range(len(ids)) if ids[i] != combo])
        test_idx = np.array([i for i in range(len(ids)) if ids[i] == combo])
        idx_tuple = (train_idx, test_idx)
        splits.append(idx_tuple)
    
    return splits


def timepoint_held_out_split(df):
    """
    Train-test split where 1 timepoint is held out in each split

    Args:
        df : Dataframe with metadata attached
    
    Returns:
        splits : List of tuples with (train_idx, test_idx)
    """

    # Store drug IDs and select unique combinations
    times = df["timepoint"].to_numpy()
    unique_times = np.unique(times)

    # Store splits
    splits = []

    for time in unique_times:
        
        # Get idx
        train_idx = np.array([i for i in range(len(times)) if times[i] != time])
        test_idx = np.array([i for i in range(len(times)) if times[i] == time])
        idx_tuple = (train_idx, test_idx)
        splits.append(idx_tuple)
    
    return splits