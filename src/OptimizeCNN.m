%% before run file,
%% set working directory to the "CNN-based_OFDM_Chaneel_Estimation"
clear;
currentDir = pwd; 
fprintf('Current Directory: %s\n', currentDir);

%% Prepare data
loadTrainData = True;
if loadTrainData
    load('train_data/trainData.mat')
else
    [trainData,trainLabels,MP] = generate_train_data(data_size,np,NRB);
    save('train_data/trainData.mat','trainData','trainLabels','MP')
end

%%
XTrain = trainData;
YTrain = trainLabels;
data_size = size(YTrain,4);
idx = randperm(data_size,data_size*0.2);
idx_test = idx(1:length(idx)/2);
idx_validation = idx(length(idx)/2+1:end);
%%
XTrain(:,:,:,idx) = [];
XTest = trainData(:,:,:,idx_test);
XValidation = trainData(:,:,:,idx_validation);

%%
YTrain(:,:,:,idx) = [];
YTest = trainLabels(:,:,:,idx_test);
YValidation = trainLabels(:,:,:,idx_validation);
