%% before run file,
%% set working directory to the "CNN-based_OFDM_Chaneel_Estimation\src"
%% run addpath("module","model","train_data")
clear;
close all;
currentDir = pwd; 
fprintf('Current Directory: %s\n', currentDir);

%% Prepare data

test_data_ratio = 0.8;
data_size=12800;  
NRB = 20;
pos=2; % pilot position pattern 1 or 2

switch pos
    case 1
        np=[3 12];
    case 2
        np=[3 8 12];
    otherwise
        warning('Unexpected pos value')
end

loadTrainData = false;
if loadTrainData
    load('..\train_data\trainData.mat')
else
    [trainData,trainLabels,MP] = generate_train_data(data_size,np,NRB);
    save('..\train_data\trainData.mat','trainData','trainLabels','MP')
end


%%
% Split real and imaginary data, then concatenate
trainData = cat(4,trainData(:,:,1,:),trainData(:,:,2,:));
trainLabels = cat(4,trainLabels(:,:,1,:),trainLabels(:,:,2,:));

batchSize=128;
% Split into training and validation sets
XValidation = trainData(:,:,:,1:batchSize);
YValidation = trainLabels(:,:,:,1:batchSize);
XTrain = trainData(:,:,:,batchSize+1:end);
YTrain = trainLabels(:,:,:,batchSize+1:end);


%% Choose Variables to Optimize
optimVars = [
    % optimizableVariable('SectionDepth',[1 3],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-4 1],'Transform','log')
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];

ObjFcn = makeObjFcn(XTrain,YTrain,XValidation,YValidation);

%% Perform Bayesian Optimization
BayesObject = bayesopt(ObjFcn,optimVars, ...
    'MaxTime',1*60*60, ...
    'IsObjectiveDeterministic',false, ...
    'UseParallel',false);

%% Evaluate the best network
% validation error
bestIdx = BayesObject.IndexOfMinimumTrace(end);
fileName = BayesObject.UserDataTrace{bestIdx};
savedStruct = load(fileName);
valError = savedStruct.valError

% test error(optional)

%% Define the objective function to optimize
function ObjFcn = makeObjFcn(XTrain,YTrain,XValidation,YValidation)
ObjFcn = @valErrorFun;
    function [valError,cons,fileName] = valErrorFun(optVars)

        % Should be specified by a global variable
        pos = 2; 
        MP = [1:2:239];
        np=[3 8 12];
        
        % Define the CNN structure

        switch pos
            case 1
                transposedCNN= transposedConv2dLayer([8 9],1,...
                    "Stride",[2 7],"Cropping",[3 1]);
            case 2
                transposedCNN= transposedConv2dLayer([8 8],1,...
                    "Stride",[2 5],"Cropping",[3 2]) ;
            otherwise
                warning('Unexpected pos value')
        end

        d=56;s=12;m=4;
        Mapping=[];
        for i=1:m
            Mapping=[
                    Mapping
                    convolution2dLayer(3,s,'Padding','same')
                    reluLayer
            ];
        end
        
        layers = [ ...
                imageInputLayer([length(MP) length(np) 1],'Normalization','none')
                convolution2dLayer(5,d,'Padding','same')
                reluLayer

                convolution2dLayer(1,s,'Padding','same')
                reluLayer
        
                Mapping
        
                convolution2dLayer(1,d,'Padding','same')
                reluLayer
        
                transposedCNN

        ];

        batchSize = 16;
        % 5 validation per epoch
        validationFrequency = round(size(XTrain,4)/batchSize /5);
        options = trainingOptions('adam', ...
            'InitialLearnRate',optVars.InitialLearnRate, ...
            'MaxEpochs',10, ...         
            'MiniBatchSize',batchSize, ...
            'L2Regularization',optVars.L2Regularization, ...
            'Shuffle','every-epoch', ...
            'Verbose',false, ...
            'Plots','training-progress', ...
            'ValidationData',{XValidation,YValidation}, ...
            'ValidationFrequency',validationFrequency, ...
            'ValidationPatience',5);

        lossFunction = "mean-squared-error";

        trainedNet =  trainnet(XTrain, YTrain, layers, lossFunction, options);
        close(findall(groot,'Tag','NNET_CNN_TRAININGPLOT_UIFIGURE'))

        % Calculate the validation error
        YPredicred = predict(trainedNet,XValidation);
        valError = YValidation - YPredicred;
        valError = mean(valError(:).^2); % need to remove dmrs_loc from error calculation
        disp(valError)

        fileName = num2str(valError) + ".mat";
        % save(fileName,'trainedNet','valError','options')
        directoryPath = '../log/test_';  
        save((directoryPath + fileName), 'trainedNet', 'valError', 'options')
        cons = [];
        
    end
end

