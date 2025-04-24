clear;
close all;
addpath("../module")
addpath("../results")
addpath("../"); % config.mのパスを追加
params = configs();

global pos np MP output_folder

%% Manually set parameters
testMSE = true;

output_folder = '2_0424';
NRB = 20;
pos=2; % pilot position pattern 1 or 2


%% Create results/output_folder directory
results_dir = fullfile('..', 'results', output_folder);
if exist(results_dir, 'dir')
    rmdir(results_dir, 's');
end
mkdir(results_dir);

% Start logging to a text file
diary(fullfile(results_dir, 'log.txt'));
diary on

disp('Parameters:')
disp(params)
disp('--------------------')

%% prepare training data
switch pos
    case 1  
        np=[3 12];
    case 2 
        np=[3 8 12];
    otherwise 
        warning('Unexpected pos value')
end
data_size=params.DataSize; 
[trainData,trainLabels,MP] = generate_train_data(data_size,np,NRB);

% Split real and imaginary data, then concatenate
trainData = cat(4,trainData(:,:,1,:),trainData(:,:,2,:));
trainLabels = cat(4,trainLabels(:,:,1,:),trainLabels(:,:,2,:));

% Split into training and validation sets
val_data_size = 2 * params.ValDataSize; %2=Real&imaginary
valData = trainData(:,:,:,1:val_data_size);
valLabels = trainLabels(:,:,:,1:val_data_size);
trainData = trainData(:,:,:,val_data_size+1:end);
trainLabels = trainLabels(:,:,:,val_data_size+1:end);


%% Choose Variables to Optimize
optimVars = [
    optimizableVariable('InitialLearnRate', [params.InitialLearnRate/100, params.InitialLearnRate*100], 'Transform', 'log')
    optimizableVariable('L2Regularization', [params.L2Regularization/100, params.L2Regularization*100], 'Transform', 'log')
    optimizableVariable('MiniBatchPower',[3 11],'Type','integer')  % 2^3=8 ～ 2^11=2048 # defalt 2^7=128
    optimizableVariable('GradientDecayFactor',[0.6 0.99])
];


disp('Optimization Variables:')
for i = 1:numel(optimVars)
    disp(optimVars(i))
end
disp('--------------------')


ObjFcn = makeObjFcn(trainData,trainLabels,valData,valLabels,params);


%% Perform Bayesian Optimization
BayesObject = bayesopt(ObjFcn,optimVars, ...
    'MaxTime',10*60, ... % seconds
    'IsObjectiveDeterministic',false, ...
    'UseParallel',false);

%% Evaluate the best network
% validation error
bestIdx = BayesObject.IndexOfMinimumTrace(end);
fileName = BayesObject.UserDataTrace{bestIdx};
matname_bestnetwork = fullfile("../results", output_folder, fileName)
% savedStruct = load(fullfile("../results", output_folder, fileName));
% valError = savedStruct.valError


%% Define the objective function to optimize
function ObjFcn = makeObjFcn(XTrain,YTrain,XValidation,YValidation, params)
ObjFcn = @valErrorFun;
    function [valError,cons,fileName] = valErrorFun(optVars)

        global pos np MP output_folder
        
        
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


        d = params.D; 
        s = params.S; 
        m = params.M;
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

        % 5 validation per epoch
        validationFrequency = round(size(XTrain,4)/params.MiniBatchSize /5);
        options = trainingOptions('adam', ...
            'MaxEpochs',params.MaxEpochs,... %params.MaxEpochs, ...
            'MiniBatchSize',2^optVars.MiniBatchPower, ...
            'InitialLearnRate',optVars.InitialLearnRate, ...
            'L2Regularization',optVars.L2Regularization, ...
            'Shuffle',params.Shuffle, ...
            'Verbose',false, ...             %'Plots','training-progress', ...
            'ValidationData',{XValidation,YValidation}, ...
            'ValidationFrequency',validationFrequency,...
            'GradientDecayFactor',optVars.GradientDecayFactor,...
            'ExecutionEnvironment','gpu');

        lossFunction = params.lossFunction;

        trainedNet =  trainnet(XTrain, YTrain, layers, lossFunction, options);
        close(findall(groot,'Tag','NNET_CNN_TRAININGPLOT_UIFIGURE'))

        % Calculate the validation error
        YPredicred = predict(trainedNet,XValidation);
        valError = YValidation - YPredicred;
        valError = mean(valError(:).^2); % need to remove dmrs_loc from error calculation

        fileName = num2str(valError) + ".mat";
        save(fullfile("../results", output_folder,fileName), 'trainedNet', 'valError', 'options')
        cons = [];
        
    end
end


%% test MSE
if testMSE
    %% fixed parameters
    numTrials=params.numTrials;  
    slots=20; % for ideal & practical LMMSE
    scs = 60e3; % subcarrier spacing 
    % MPと一緒ならmp = MPと書き換え
    mp = [1:2:12]; % pilot carrier of 1RB
    Nfft=512;
    % M=4;bps=log2(M);avp=0.5; % bpsk
    M=16;bps=log2(M);avp=2.5; % QPSK
    Nslot = 1; % number of RBs along time axis
    CP = 0.07; %　percentage of cyclic prefic
    m_1RB = 12; n_1RB =14; % 1RB= m_1RB carriers * n_1RB time slots
    fc=50*10^9; % carrier freqency
    
    % if pos==1
    %     np=[3 12]; % pilot timeslot of 1RB
    % elseif pos==2
    %     np=[3 8 12];
    % end
    
    m_1user= m_1RB*NRB; n_1user=n_1RB*Nslot;
    Ndmrs=NRB*length(mp)* Nslot*length(np); 
    Ns=m_1user*n_1user; % data+dmrs symbol
    % Number of samples transmitted by base station
    Ng=ceil(Nfft*0.07);Nofdm=Nfft+Ng; % samples/ofdm symbol
    Tofdm= Nofdm/(scs*Nfft);% OFDM symbol duration
    
    
    %% dmrs location
    MP=[];NP=[];
    for i=1:NRB MP=[MP,(i-1)*m_1RB+mp]; end
    for i=1:Nslot NP=[NP,(i-1)*n_1RB+np]; end
    a = repmat((NP-1)*m_1user,length(MP),1);
    b = repmat(MP',1,length(NP));
    c = a+b;
    dmrs_loc = c(:);
    
    
    %%  DMRS genaration
    dmrsBit=randi([0 1],1,bps*(Ndmrs/length(NP)));
    dmrsSym_m=modu(dmrsBit,M);
    dmrsSym=[]; 
    for i=1:length(NP) dmrsSym=[dmrsSym, dmrsSym_m]; end

    %% Monte preparing
    SNRdB = 0:5:30;
    for k=1:length(SNRdB)
        if  SNRdB(k)>=25; Monte=15*numTrials;
        elseif SNRdB(k)>=15; Monte=10*numTrials;
        elseif SNRdB(k)>=10; Monte=6*numTrials;
        else Monte=3*numTrials;
        end;
    
    
        fprintf('%d[dB] %d Simulations\n',SNRdB(k),Monte);
    
        for j=1:Monte
            %% randomness
            t0=randi([0 10000],1,1); % randam initail time (1=1ofdm symbol duration)
            v=randi([0 60],1,1);
            DelayProfile = char(randsample(["TDL-A","TDL-B","TDL-C",...
                ], 1));
            DS=randsample([55 228],1,1);
            fd = v*1000/3600 * fc/(3*10^8); 
    
            %% RG
            dataBit = randi([0 1],1,bps*(Ns-Ndmrs));dataSym = modu(dataBit,M);
            % dmrsSym=10*[1:Ndmrs];dataSym=[1:Ns-Ndmrs];
            [RG,RGseq]=GenRG(Nfft,m_1RB,n_1RB, ...
                NRB,Nslot,...
                dataSym,dmrsSym,dmrs_loc);
            
            %% IFFT,add GI
            IFFTout = ifft(RG);GI=IFFTout((end-Ng+1):end,:);
            sOFDM=([GI; IFFTout]);%sOFDM=sOFDM(:); % vectorization
    
            %% Channel 
            h = nrtdl(fd,Tofdm,n_1user,t0,DelayProfile,DS);
            Ndelay=size(h,1);
            chOut = zeros(length(sOFDM(:)) + Ndelay-1 ,1);
            for i = 1:n_1user
                S = (i-1)*size(sOFDM,1)+1; E = S+ size(sOFDM,1)-1 + Ndelay-1;
                chOut(S:E) =chOut(S:E)+ conv(sOFDM(:,i),h(:,i)); % Channel path (convolution)
            end
            chOut(length(sOFDM(:)) +1:end)=[];
            awgn = randn(size(chOut)) + 1i*randn(size(chOut)); 
            Pn = 10^(-SNRdB(k)/10)/Nfft *avp;%/Rofdm
            chOut = chOut + awgn*sqrt(Pn/2);
            
            %% remove GI, FFT
            rxSPCOut = reshape(chOut,Nofdm,n_1user);
            rxSPCOut (1:Ng,:) = []; 
            FFTout = fft(rxSPCOut);
            FFTout = FFTout(1:m_1user,:);
            rRGseq = transpose(FFTout(:));
            
            
            %% estimation
            H_perfect=fft(h,Nfft,1); %列(遅延)方向にfft
            H_perfect=H_perfect(1:m_1user,:);
    
    
            %% LS
            H_LS = rRGseq(dmrs_loc)./dmrsSym; 
            H_LS = reshape(H_LS,length(MP),length(NP));    
     
            %% FSRCNN(bays)
            load(matname_bestnetwork);
            nnInput = cat(4,real(H_LS),imag(H_LS));
            H_cnn = predict(trainedNet,nnInput);
            H_FSRCNN_bays(1:m_1user,1:n_1user) = H_cnn(:,:,1,1)+ 1i*H_cnn(:,:,1,2);
            H_e = H_perfect-H_FSRCNN_bays;
            H_e=H_e(:); H_e(dmrs_loc)=[];
            mse(j) = mean(abs(H_e).^2)  ;
    
        end
        MSE_FSRCNN_bays(k)=mean(mse,2);
    end     
    
    save(fullfile(results_dir, "MSE_bays.mat"),"MSE_FSRCNN_bays")
    %% plot MSE
    figure;
    markersize=15;
    load("../results/MSE_traditional.mat");
    load("../results/MSE_FSRCNN_normal.mat");
    semilogy(SNRdB,MSE(1,:),'r+-','MarkerSize',markersize);hold on;
    semilogy(SNRdB,MSE(2,:),'go-','MarkerSize',markersize);hold on;grid on;
    semilogy(SNRdB,MSE(3,:),'bo-','MarkerSize',markersize);hold on;
    semilogy(SNRdB,MSE_FSRCNN_normal,'ksquare-','MarkerSize',markersize);hold on;
    semilogy(SNRdB,MSE_FSRCNN_bays,'k^-','MarkerSize',markersize);

    xlabel('SNR[dB]') 
    ylabel('MSE')
    legend('LS',...
        'ideal LMMSE','practical LMMSE',...
        '深層学習（通常）','深層学習（bays）',...
        'FontSize',22);
    xlim([0 SNRdB(end)])
    set(gca,'FontSize',22)
    
end

toc;
    
     addpath("C:\Users\ide-lab\OneDrive - Kyushu University\研究2025\code\CNN-based_OFDM_Chaneel_Estimation\results")
    





