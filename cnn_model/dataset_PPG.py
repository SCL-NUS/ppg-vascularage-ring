import pandas as pd
import numpy as np
from torch.utils.data import Dataset
import torch
from torch.utils.data import DataLoader



class PPGDataset(Dataset):
    def __init__(self, file_path):
        data = pd.read_csv(file_path, header=None)
        # Assuming the first row contains metadata and should be excluded from the dataset
        self.X = torch.tensor(data.iloc[:, :-2].values, dtype=torch.float32)  # Excludes the last two columns and the first row
        self.y = torch.tensor(data.iloc[:, -2].values, dtype=torch.float32)  # Takes the second last column as labels

    def __len__(self):
        # Length is determined by the size of X or y, both should be the same
        return len(self.X)

    def __getitem__(self, idx):
        # Provide access to a single item pair (feature, label)
        if idx >= len(self.X):
            raise IndexError("Index out of bounds")  # This will explicitly prevent accessing an invalid index
        return self.X[idx], self.y[idx]

