import torch
from torch.utils.data import DataLoader
import torch.optim as optim
import torch.nn as nn
from model import AgeEstimationCNN
from dataset_PPG import PPGDataset
import csv
import numpy as np

def train_model(model, train_loader, val_loader, criterion, optimizer, epochs, patience):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(device)
    model = model.to(device)
    best_val_loss = 1e8
    epochs_no_improve = 0
    losses = {'train': [], 'val': []}
    train_rmse = []  # List to save RMSE for each training epoch
    val_rmse = []  # List to save RMSE for each validation epoch
    train_mae = []  # List to save MAE for each training epoch
    val_mae = []  # List to save MAE for each validation epoch
    train_r2 = []  # List to save R² for each training epoch
    val_r2 = []  # List to save R² for each validation epoch

    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        all_pred_train = []
        all_labels_train = []

        for inputs, labels in train_loader:
            inputs, labels = inputs.to(device), labels.to(device)

            optimizer.zero_grad()
            outputs = model(inputs.unsqueeze(1))
            loss = criterion(outputs.squeeze(), labels)
            loss.backward()
            optimizer.step()

            running_loss += loss.item()
            all_pred_train.append(outputs.detach().view(-1).cpu().numpy())
            all_labels_train.append(labels.detach().cpu().numpy())

        train_loss = running_loss / len(train_loader)
        losses['train'].append(train_loss)
        all_pred_train = np.concatenate(all_pred_train)
        all_labels_train = np.concatenate(all_labels_train)

        # all_pred_train = torch.tensor(all_pred_train)
        # all_labels_train = torch.tensor(all_labels_train)
        rmse_train = np.sqrt(np.mean((all_pred_train - all_labels_train) ** 2))
        train_rmse.append(rmse_train.item())

        mae_train = np.mean(np.abs(all_pred_train - all_labels_train))
        train_mae.append(mae_train.item())

        r2_train = 1 - np.sum((all_pred_train - all_labels_train) ** 2) / np.sum((all_labels_train - np.mean(all_labels_train)) ** 2)
        train_r2.append(r2_train.item())

        all_pred_val = []
        all_labels_val = []
        with torch.no_grad():
            model.eval()
            val_loss = 0
            for inputs, labels in val_loader:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs.unsqueeze(1))
                loss = criterion(outputs.squeeze(), labels)
                val_loss += loss.item()

                all_pred_val.append(outputs.detach().view(-1).cpu().numpy())
                all_labels_val.append(labels.detach().cpu().numpy())

            val_loss /= len(val_loader)
            losses['val'].append(val_loss)

            # all_pred_val = torch.tensor(all_pred_val)
            # all_labels_val = torch.tensor(all_labels_val)
            # rmse_val = torch.sqrt(torch.mean((all_pred_val.float() - all_labels_val.float()) ** 2))
            # val_rmse.append(rmse_val.item())
            all_pred_val = np.concatenate(all_pred_val)
            all_labels_val = np.concatenate(all_labels_val)
            rmse_val = np.sqrt(np.mean((all_pred_val - all_labels_val) ** 2))
            val_rmse.append(rmse_val.item())
            mae_val = np.mean(np.abs(all_pred_val - all_labels_val))
            val_mae.append(mae_val.item())

            r2_val = 1 - np.sum((all_pred_val - all_labels_val) ** 2) / np.sum((all_labels_val - np.mean(all_labels_val)) ** 2)
            val_r2.append(r2_val.item())

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                epochs_no_improve = 0
            else:
                epochs_no_improve += 1
            if epochs_no_improve == patience:
                break

    return losses, train_rmse, val_rmse, train_mae, val_mae, train_r2, val_r2

def save_results_to_csv(losses, train_rmse, val_rmse, train_mae, val_mae, train_r2, val_r2, train_file, val_file):
    with open(train_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Epoch', 'Train Loss', 'Train RMSE', 'Train MAE', 'Train R2'])
        for epoch, (loss, rmse, mae, r2) in enumerate(zip(losses['train'], train_rmse, train_mae, train_r2)):
            writer.writerow([epoch, loss, rmse, mae, r2])

    with open(val_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Epoch', 'Validation Loss', 'Validation RMSE', 'Validation MAE', 'Validation R2'])
        for epoch, (loss, rmse, mae, r2) in enumerate(zip(losses['val'], val_rmse, val_mae, val_r2)):
            writer.writerow([epoch, loss, rmse, mae, r2])