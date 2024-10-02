clear
trainModel =true;

NRB = 20; % 1RB use 12 subcarrier 
pos=1;

matName_FSRCNN=sprintf("FSRCNN_Pos%d_RB%d",pos,NRB);
switch pos
    case 1
        np=[3 12];
        transposedCNN= transposedConv2dLayer([8 9],1,...
             "Stride",[2 7],"Cropping",[3 1]);
    case 2
        np=[3 8 12];
        transposedCNN= transposedConv2dLayer([8 8],1,...
             "Stride",[2 5],"Cropping",[3 2]) ;
     otherwise
        warning('Unexpected pos value')
end

if trainModel
    % Generate the training data
    [trainData,trainLabels,MP] = TrainData_FSRCNN(12800,np,NRB);

    % Set the number of examples per mini-batch
    batchSize=128;

    % Split real and imaginary grids into 2 image sets, then concatenate
    trainData = cat(4,trainData(:,:,1,:),trainData(:,:,2,:));
    trainLabels = cat(4,trainLabels(:,:,1,:),trainLabels(:,:,2,:));
    
    % Split into training and validation sets
    valData = trainData(:,:,:,1:batchSize);
    valLabels = trainLabels(:,:,:,1:batchSize);
    trainData = trainData(:,:,:,batchSize+1:end);
    trainLabels = trainLabels(:,:,:,batchSize+1:end);
    
    % Validate roughly 5 times every epoch
    % iteration15/5 =3 ⇨1回/3iteration,5回/epoch
    valFrequency = round(size(trainData,4)/batchSize/5);

    % Define the CNN structure
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
    
    % Set up a training policy
    options = trainingOptions('adam', ...
        'InitialLearnRate',1e-3, ...
        'MaxEpochs',50, ...
        'Shuffle','every-epoch', ...
        'Verbose',false, ...
        'Plots','training-progress', ...
        'MiniBatchSize',batchSize, ...
        'ValidationData',{valData, valLabels}, ...
        'ValidationFrequency',valFrequency, ...
        'ValidationPatience',5);
    
    % Train the network. The saved structure trainingInfo contains the
    % training progress for later inspection. This structure is useful for
    % comparing optimal convergence speeds of different optimization
    % methods.
    [channelEstimationCNN,trainingInfo] = trainNetwork(trainData, ...
        trainLabels,layers,options);


    save(matName_FSRCNN, 'channelEstimationCNN')
end


