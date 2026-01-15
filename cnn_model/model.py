import torch
import torch.nn as nn
import torch.nn.functional as F

class AgeEstimationCNN(nn.Module):
    def __init__(self):
        super(AgeEstimationCNN, self).__init__()
        # self.conv1 = nn.Conv1d(in_channels=1, out_channels=16, kernel_size=10, stride=1, padding=5)
        # self.conv2 = nn.Conv1d(in_channels=16, out_channels=32, kernel_size=8, stride=1, padding=4)
        self.conv1 = nn.Conv1d(in_channels=1, out_channels=16, kernel_size=11, stride=1, padding=5)
        self.conv2 = nn.Conv1d(in_channels=16, out_channels=32, kernel_size=9, stride=1, padding=4)
        # self.fc1 = nn.Linear(32 * 184, 1024)
        self.fc1 = nn.Linear(32 * 200, 1024)
        self.fc2 = nn.Linear(1024, 1024)
        self.out = nn.Linear(1024, 1)

        self.dropout = nn.Dropout(0.2)

    def forward(self, x):
        x = F.relu(self.conv1(x))# Applying the first convolutional layer followed by ReLU activation
        x = F.relu(self.conv2(x)) # Applying the second convolutional layer followed by ReLU activation
        x = x.view(x.size(0), -1) # Flatten the output for the fully connected layer
        x = self.dropout(F.relu(self.fc1(x)))# First fully connected layer with ReLU activation and dropout
        x = self.dropout(F.relu(self.fc2(x)))# Second fully connected layer with ReLU activation and dropout
        x = self.out(x)
        return x
