# Functions for splitting train/test set according to specified metrics
import numpy as np
import pandas as pd
import sklearn
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

def condition_to_replicate_idx(condition_list):
    """
    Function to convert a list of conditions into replicate labels

    Args:
        condition_list : List of condition labels, ex: ["12CIP1hr-a", "12CIP1hr-b",...]

    Returns:
        idx_list : List of condition labels converted into replicate indices, ex: [0, 0, 1,...]
    """
    idx_list = [0]

    for i in range(1, len(condition_list)):
        prev_condition = condition_list[i-1]   
        curr_condition = condition_list[i]

        # Exclude last character since "-a" "-b"
        if prev_condition[:-1] == curr_condition[:-1]:
            idx_list.append(idx_list[i-1])
        
        else:
            idx_list.append(idx_list[i-1] + 1)
    
    return idx_list
    
def find_first_alpha(str):
    """
    Function to find the index of the the first letter in a string
    (will be used for drug name parser function)

    Args:
        str : string

    Returns:
        first_alpha : idx of first letter
    """
    for i in range(len(str)):
        if str[i].isalpha():
            return i
            break 

def condition_to_drug(condition_label):
    """
    Function to convert a condition label to a drug name

    Args:
        condition_label : Condition label, ex: "12CIP1hr-a"

    Returns:
        drug_name : Drug name, ex: "CIP"
    """
    # Single drug case
    if len(condition_label) < 11:

        # If NDC sample
        if "NDC" in condition_label:
            return "NDC"
        
        else:

            # Find drug name using letter search
            first_alpha_idx = find_first_alpha(condition_label)
            drug_name = condition_label[first_alpha_idx:first_alpha_idx + 3]
            return drug_name
    
    # Multiple drug case
    else:

        # Find first letter and extract drug name
        first_idx = find_first_alpha(condition_label)
        drug1 = condition_label[first_idx:first_idx + 3]

        # Find first letter in remaining string and extract drug name
        substr     = condition_label[first_idx + 3:]
        second_idx = find_first_alpha(substr)
        drug2      = substr[second_idx:second_idx + 3]

        # Merge drug names
        drug_name = drug1 + "+" + drug2

        return drug_name
    
def custom_train_test_split(data_df, test_size):
    """
    Function to generate custom train test split according to 2 rules
        1. Keep triplicates in the same set
        2. Maintain equal drug class distribution across train/test

    Args:
        data_df     [N,G+1] : Dataframe with all TPM values (N samples, G genes + 1 column of CFUs)
                              Index should contain condition names
        test_size           : 0.1, 0.2, 0.25, or 0.33 for test set percentage

    Returns:
        X_train, y_train, X_test, y_test
    """
    X = data_df.iloc[:, data_df.columns != "CFU"]
    y = data_df.iloc[:, data_df.columns == "CFU"]

    condition_names = list(X.index)

    # Class and replicate labels
    drug_labels = np.array([condition_to_drug(cond) for cond in condition_names])
    rep_labels  = np.array(condition_to_replicate_idx(condition_names))

    # Convert test_size to arg for groupkfold
    n_splits = round(1 / test_size)

    # Split, preserving replicates and class distribution
    splitter = StratifiedGroupKFold(n_splits = n_splits, shuffle = True, random_state = 112)
    train_idx, test_idx = next(splitter.split(X, y = drug_labels, groups = rep_labels)) # Take first set of idx
    X_train, y_train = X.iloc[train_idx], y.iloc[train_idx]
    X_test, y_test   = X.iloc[test_idx], y.iloc[test_idx]

    return X_train, y_train, X_test, y_test

def apply_scaler(X_train, X_test):
    """
    Function to train scaler on training set, then apply to test

    Args:
        X_train, X_test

    Returns:
        X_scaled_train, X_scaled_test
    """
    scaler = StandardScaler()
    scaler.set_output(transform = "pandas")

    scaler.fit(X_train)

    X_scaled_train = scaler.transform(X_train)
    X_scaled_test  = scaler.transform(X_test)
    
    return X_scaled_train, X_scaled_test
 
def full_train_test_split(data_df, test_size):
    """
    Function to run full train-test split procedure, including grouped split and scaler

    Args:
        data_df     [N,G+1] : Dataframe with all TPM values (N samples, G genes + 1 column of CFUs)
                              Index should contain condition names
        test_size           : 0.1, 0.2, 0.25, or 0.33 for test set percentage

    Returns:
        X_train, y_train, X_test, y_test : Standard scaled
    """
    X_train_unscaled, y_train, X_test_unscaled, y_test = custom_train_test_split(data_df = data_df, test_size = test_size)
    X_train, X_test = apply_scaler(X_train_unscaled, X_test_unscaled)

    return X_train, y_train, X_test, y_test