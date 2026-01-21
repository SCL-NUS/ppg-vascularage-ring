import torch
import matplotlib.pyplot as plt
import numpy as np
from torchcam.methods import GradCAM
from torch.utils.data import DataLoader
from model import AgeEstimationCNN
from dataset_PPG_withID import PPGDataset
import os
import pandas as pd
import re


def parse_model_string(model_string):
    # Adjust the regex to match the given pattern more accurately
    match = re.match(r"(\w+_\w+)_model_(bs\d+_lr\d+\.\d+_ep\d+)_patience(\d+)", model_string)
    if match:
        model_type = match.group(1)
        parameters = f"{match.group(2)}_patience{match.group(3)}"
        # fold_number = match.group(4)
        return model_type, parameters
    return None, None




def apply_grad_cam(model, input_tensor, target_layer):
    # Initialize the Grad-CAM object for the specified layer
    model.eval()
    cam_extractor = GradCAM(model, target_layer=target_layer)
    # print(model)
        # Ensure the model is in evaluation mode
   
    
    # Forward pass through the model to get the output
    output = model(input_tensor)

    
    # Generate the Grad-CAM mask for the regression output
    activation_map = cam_extractor(0, output)  # Targeting the entire output, assuming it's a scalar
 
    
 # Convert the activation map to numpy
    activation_map = activation_map[0].cpu().detach().numpy()

    # Flatten the activation map (assuming it's not already flat)
    activation_map = activation_map.flatten()
    
    # # Interpolate the activation map to match the length of the PPG pulse
    # target_length = input_tensor.shape[-1]  # Get the length of the PPG pulse from the input tensor
    # activation_map = np.interp(
    #     np.arange(target_length),  # Target length (e.g., 200 for your PPG pulse)
    #     np.linspace(0, target_length - 1, num=activation_map.shape[0]),  # Activation map indices
    #     activation_map  # Activation map values
    # )
    output = output.detach().cpu().numpy()


    return activation_map, output



def average_heatmaps_per_subject(test_dataset, model, target_layer='conv2'):

    # Initialize an empty DataFrame with appropriate columns
    columns = ['subject_id'] + [f'heatmap_{i}' for i in range(200)] + [f'ppg_{i}' for i in range(200)] + ['actual_age', 'predicted_age']
    df = pd.DataFrame(columns=columns)
# Loop through the dataset
    for i in range(len(test_dataset)):
        ppg_sample, label, subject_id = test_dataset[i][0], test_dataset[i][1], test_dataset[i][2]
        ppg_sample = ppg_sample.unsqueeze(0).unsqueeze(0)  # Shape: [1, 1, 200]
        subject_id=subject_id.squeeze().cpu().numpy()
        label=label.squeeze().cpu().numpy()
        # Apply Grad-CAM
        target_layer='conv2'
        activation_map, predict = apply_grad_cam(model, ppg_sample, target_layer)

        # Prepare the data for the DataFrame
        heatmap_row = activation_map.flatten()  # Flatten heatmap to 1D (if it's not already)
        ppg_row = ppg_sample.squeeze().cpu().numpy().flatten()  # Flatten PPG to 1D
        row = [subject_id.item()] + heatmap_row.tolist() + ppg_row.tolist() + [label.item(), predict.item()]

        # Append the data to the DataFrame
        df.loc[len(df)] = row
    return df


############################################################################################################
############################################################################################################

curr_dir='/directory where scripts are' 
data_main_dir='/directory where individual pulse waveforms are saved' 
DS_pth =os.path.join(data_main_dir, 'FT')
data_dir = os.path.join(DS_pth , 'folds')



output_dir = os.path.join(curr_dir,'GradCam','grad_cam_visualizations_avg_allfolds')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


df = pd.read_csv(os.path.join(curr_dir,'best_hyperparameters_FT.csv')) # read table with models with best hyperparameters 

md_src = '/directory where trainedl models are stored' # example name string of a trained model: trained_model_bs128_lr0.001_ep200_patience100_fold_1.pth


complete_df = pd.DataFrame()


for fold_number in range(10):

    model_string = df['Model_Name']
    mdl_str=str(model_string[0])
    model_type, parameters = parse_model_string(mdl_str)
    md=model_type 
    md_srcpath = os.path.join(md_src, md)
    # criterion = nn.MSELoss() if 'MSE' in md else nn.L1Loss()


    test_file = os.path.join(data_dir, f'fold{int(fold_number+1)}', 'test.csv') # get waveforms used in test set
    test_dataset = PPGDataset(test_file)
    batch_size = len(test_dataset)
    # batch_size = 20
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False, drop_last=False)
    model = AgeEstimationCNN()
    dataloader = test_loader


    model_file=f'trained_model_{parameters}_fold_{fold_number+1}.pth'
    model_path = os.path.join(md_srcpath, model_file)
    print(model_path)
    model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))

    # Set the model to evaluation mode
    model.eval()  # Place this line here
    ds= average_heatmaps_per_subject(test_dataset, model)
    grouped_df = ds.groupby('subject_id').mean().reset_index()
    complete_df = pd.concat([complete_df, grouped_df], ignore_index=True)

    print("Average heatmaps per subject")

# Save the DataFrame to a file if needed
complete_df.to_csv(os.path.join(curr_dir,'GradCam','Allsubject_data2.csv'), index=False)    


# Plot

for index, row in complete_df.iterrows():
        subject_id = row['subject_id']
        average_pulse = row[[f'ppg_{i}' for i in range(200)]].values  # Extract averaged PPG pulse
        average_age = row['actual_age']
        average_pred = row['predicted_age']
        heatmap = row[[f'heatmap_{i}' for i in range(200)]].values  # Extract averaged heatmap

        # Plot the averaged PPG pulse
        plt.figure(figsize=(10, 4))
        plt.plot(average_pulse, label='Average PPG Signal', color='blue')
        
        # Plot the averaged heatmap
        plt.imshow(heatmap[np.newaxis, :], cmap='jet', aspect='auto', alpha=0.5,
                extent=[0, len(average_pulse), average_pulse.min(), average_pulse.max()])

        # Add colorbar and title
        plt.colorbar(label='Activation')
        plt.title(f'Subject {subject_id}, Actual Age: {average_age:.2f}, Predicted Age: {average_pred:.2f}')
        print(f'Subject {subject_id}')
        
        # Save the figure
        plt.savefig(os.path.join(output_dir, f'grad_cam_average_subject_{subject_id}_Age_{average_age}.png'))
        plt.close()



print("done")

