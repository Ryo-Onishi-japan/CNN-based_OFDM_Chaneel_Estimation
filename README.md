## Implementation of the paper ["An Evaluation of Channel Estimation Method Using Deep Learning for OFDM System"](https://ieeexplore.ieee.org/document/10625125)  
Presented in: 2024 Fifteenth International Conference on Ubiquitous and Future Networks (ICUFN)

## Abstract

In the time- and frequency-variant mobile radio channel such as the Fifth- and Fourth-Generation mobile communications systems (5G, 4G), it is very important to estimate the channel coefficients (gains). We propose an estimation method of channel coefficients using deep learning. The high accuracy of estimation can be evaluated using the deep learning method for super-resolution (SR) Network compared with the conventional method. The estimation of channel gains using orthogonal frequency-division multiplexing (OFDM) can be replaced with the SR problem and applied to the SR Network method. We show that the accuracy of the proposed method is higher than the conventional methods.

## Presentation slides
- [short (Japanese)](https://1drv.ms/p/c/91b183a988eac49f/EUXK5924eRtKjSmBextnGSkBwHFdLO41UaTgX4DAh12VNg?e=8cuDKf)
- [long (English)](https://1drv.ms/p/c/91b183a988eac49f/EVNWYS9BwAxEmv601Qu1fAEBDTzdUJJrPkj1lcu5cI0g6w?e=lX507H)


## How to Run
1. Run `train_cnn.m` to train the model.
2. Run `main.m` to test the deep learning method and conventional methods.


## Major Files and Folders
- `generate_train_data.m`: Script to generate training data.
- `module/`: containing utility functions and modules used in the project.
- `model/`:containing the deep learning model definitions.
- `train_data/`: generated training data is stored.


## Notes
- MATLAB Deep Learning Toolbox is required.
