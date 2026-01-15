import argparse
import torch
import os
import pandas as pd
from training_utils4 import train_model, save_results_to_csv
import torch.optim as optim
import torch.nn as nn
from torch.utils.data import DataLoader
from model import AgeEstimationCNN
from dataset_PPG import PPGDataset
# from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts

def main():
    # Setting up ArgumentParser
    parser = argparse.ArgumentParser(description='Train a model for age estimation using cross-validation')
    parser.add_argument('--data_dir', type=str, required=True, help='Directory containing the dataset files')
    parser.add_argument('--result_path', type=str, required=True, help='Path to save the results')
    parser.add_argument('--model_path', type=str, required=True, help='Path to save the trained models')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate')
    parser.add_argument('--batch-size', type=int, default=20, help='Batch size')
    parser.add_argument('--epochs', type=int, default=900, help='Number of epochs')
    parser.add_argument('--patience', type=int, default=90, help='Patience for early stopping')
    parser.add_argument('--optimizer_name', type=str, default='Adam', help='Optimizer Type')
    parser.add_argument('--loss_function', type=str, default='MSE', help='Loss Function type')
    args = parser.parse_args()


    for idx in range(10):
     

        # test_file = os.path.join(args.data_dir, f'fold{idx + 1}', 'test.csv')
        val_file = os.path.join(args.data_dir, f'fold{idx + 1}', 'validation.csv')
        train_file = os.path.join(args.data_dir, f'fold{idx + 1}', 'train.csv')

        train_dataset = PPGDataset(train_file) 
        val_dataset = PPGDataset(val_file)
        # test_dataset = PPGDataset(test_file)


################################################################

        # Concatenate training datasets
       
        train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True,drop_last=True)
        val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False,drop_last=True)

        # Model, Optimizer, Loss function setup
        model = AgeEstimationCNN()

            # Select optimizer
        if args.optimizer_name == 'Adam':
            optimizer = optim.Adam(model.parameters(), lr=args.lr,weight_decay=1e-7)
        else:
            optimizer = optim.SGD(model.parameters(), lr=args.lr, momentum=0.9,weight_decay=1e-7)


        # T_0=20
        # T_mult=2
        # scheduler = CosineAnnealingWarmRestarts(optimizer, T_0, T_mult)

        # if args.optimizer_name == 'Adam':
        if args.loss_function == 'L1':
                 criterion = nn.L1Loss()
        elif args.loss_function == 'MSE':
                criterion = nn.MSELoss()
        # elif args.optimizer_name == 'SGD':
        #     criterion = nn.L1Loss()  # Only uses L1Loss as per your specification

        # optimizer = optim.Adam(model.parameters(), lr=args.lr)
        # criterion = nn.MSELoss()

      # Training the model
        losses, train_rmse, val_rmse, train_mae, val_mae, train_r2, val_r2 = train_model(model, train_loader, val_loader, criterion, optimizer, args.epochs, args.patience)
############################################
        fold_name=f'fold_{idx + 1}'
        folder_name=f'{args.optimizer_name}_{args.loss_function}'

        model_folder = os.path.join(args.model_path, folder_name)
        if not os.path.exists(model_folder):
            os.makedirs(model_folder)
        
        # filenames for loss values that include parameter values
        # filename_suffix = f"bs{args.batch_size}_lr{args.lr}_ep{args.epochs}_{fold_name}"
        filename_suffix = f"bs{args.batch_size}_lr{args.lr}_ep{args.epochs}_patience{args.patience}_{fold_name}"
    
        # train_results_filename = os.path.join(args.result_path, f'SGD_train_loss_results_{filename_suffix}_14052023.csv')
        # val_results_filename = os.path.join(args.result_path, f'SGD_val_loss_results_{filename_suffix}_14052023.csv')
        result_folder=os.path.join(args.result_path,folder_name)
        if not os.path.exists(result_folder):
            os.makedirs(result_folder)

        
        train_results_filename = os.path.join(result_folder,f'train_loss_results_{filename_suffix}.csv')
        val_results_filename = os.path.join(result_folder,f'val_loss_results_{filename_suffix}.csv')

        model_filename = os.path.join(model_folder, f'trained_model_{filename_suffix}.pth')


        save_results_to_csv(losses, train_rmse, val_rmse, train_mae, val_mae, train_r2, val_r2, train_results_filename, val_results_filename)
        torch.save(model.state_dict(), model_filename)

if __name__ == "__main__":
    main()
