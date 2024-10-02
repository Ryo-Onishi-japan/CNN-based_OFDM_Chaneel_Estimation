## Implementation of the paper ["An Evaluation of Channel Estimation Method Using Deep Learning for OFDM System"](https://ieeexplore.ieee.org/document/10625125)  
Presented in: 2024 Fifteenth International Conference on Ubiquitous and Future Networks (ICUFN)

## Abstract

In the time- and frequency-variant mobile radio channel such as the Fifth- and Fourth-Generation mobile communications systems (5G, 4G), it is very important to estimate the channel coefficients (gains). We propose an estimation method of channel coefficients using deep learning. The high accuracy of estimation can be evaluated using the deep learning method for super-resolution (SR) Network compared with the conventional method. The estimation of channel gains using orthogonal frequency-division multiplexing (OFDM) can be replaced with the SR problem and applied to the SR Network method. We show that the accuracy of the proposed method is higher than the conventional methods.

## Presentation slides
- [short (Japanese)](https://1drv.ms/p/c/91b183a988eac49f/EUXK5924eRtKjSmBextnGSkBwHFdLO41UaTgX4DAh12VNg?e=8cuDKf)
- [long (English)](https://1drv.ms/p/c/91b183a988eac49f/EVNWYS9BwAxEmv601Qu1fAEBDTzdUJJrPkj1lcu5cI0g6w?e=lX507H)

## set-up
    cd CNN-based_OFDM_Chaneel_Estimation
    addpath("module","model","train_data")
## Data Preparation
1. Run `generate_train_data.m` to generate training data.
    ```matlab
   dataSize = 12800;
   np = [3, 8, 12]; % Pilot positions along the time axis
   NRB = 20; % Number of resource blocks along the frequency axis
   [trainData, trainLabels, MP] = generate_train_data(dataSize, np, NRB);
   % MP - Pilot positions along the frequency axis
    ```
2. Save the generated data to `train_data.mat`.
    ```matlab
    save('train_data.mat', 'train_data', 'train_labels', 'mp');
    ```

## Training
1. Run `train_cnn.m` to train the model.
    ```matlab
    [channel_estimation_cnn, training_info] = trainNetwork(train_data, train_labels, layers, options);
    save('channel_estimation_cnn.mat', 'channel_estimation_cnn');
    ```

## Model Structure
1. The model structure is defined as follows:
    ```matlab
    layers = [ ...
        imageInputLayer([length(mp) length(np) 1],'Normalization','none')
        convolution2dLayer(5, d, 'Padding', 'same')
        reluLayer
        convolution2dLayer(1, s, 'Padding', 'same')
        reluLayer
        mappingLayer
        convolution2dLayer(1, d, 'Padding', 'same')
        reluLayer
        transposedCNNLayer
        regressionLayer
    ];
    ```

## File Structure
- `generate_train_data.m`: Script to generate training data.
- `train_cnn.m`: Script to train the model.
- `optimize_cnn.m`: Script to save training data.
- `train_data.mat`: Generated training data.

## How to Run
1. Run `generate_train_data.m` to prepare the data.
2. Run `train_cnn.m` to train the model.

## Notes
- MATLAB Deep Learning Toolbox is required.
