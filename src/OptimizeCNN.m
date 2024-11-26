%% before run file,
%% set working directory to the "CNN-based_OFDM_Chaneel_Estimation"
clear;
currentDir = pwd; 
fprintf('Current Directory: %s\n', currentDir);

%% Prepare data
data_size=12800;
NRB = 20;
global pos=2; % pilot position pattern 1 or 2

switch pos
    case 1
        np=[3 12];
    case 2
        np=[3 8 12];
    otherwise
        warning('Unexpected pos value')
end

loadTrainData = true;
if loadTrainData
    load('train_data/trainData.mat')
else
    [trainData,trainLabels,MP] = generate_train_data(data_size,np,NRB);
    save('train_data/trainData.mat','trainData','trainLabels','MP')
end


%%
% Split real and imaginary data, then concatenate
trainData = cat(4,trainData(:,:,1,:),trainData(:,:,2,:));
trainLabels = cat(4,trainLabels(:,:,1,:),trainLabels(:,:,2,:));

% Split into training, validation, test sets
XTrain = trainData;
YTrain = trainLabels;
data_size = size(YTrain,4);
idx = randperm(data_size,data_size*0.2);
idx_test = idx(1:length(idx)/2);
idx_validation = idx(length(idx)/2+1:end);

XTrain(:,:,:,idx) = [];
XTest = trainData(:,:,:,idx_test);
XValidation = trainData(:,:,:,idx_validation);

YTrain(:,:,:,idx) = [];
YTest = trainLabels(:,:,:,idx_test);
YValidation = trainLabels(:,:,:,idx_validation);

%% Choose Variables to Optimize
optimVars = [
    % optimizableVariable('SectionDepth',[1 3],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-2 1],'Transform','log')
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];

ObjFcn = makeObjFcn(XTrain,YTrain,XValidation,YValidation);

%% Perform Bayesian Optimization
BayesObject = bayesopt(ObjFcn,optimVars, ...
    'MaxTime',14*60*60, ...
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
        
        % Define the CNN structure
        switch global pos
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

                regressionLayer
        ];

        batchSize = 128;
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
            'ValidationFrequency',validationFrequency);
            'ValidationPatience',5);

        lossFunction = "mean-squared-error";

        trainedNet =  trainnet(trainData, trainLabels, layers, lossFunction, options);
        close(findall(groot,'Tag','NNET_CNN_TRAININGPLOT_UIFIGURE'))

        % Calculate the validation error
        YPredicred = predict(trainedNet,XValidation);
        valError = YValidation - YPredicred;
        valError = mean(valError(:).^2); % need to remove dmrs_loc from error calculation

        fileName = num2str(valError) + ".mat";
        save(fileName,'trainedNet','valError','options')
        cons = [];
        
    end
end

