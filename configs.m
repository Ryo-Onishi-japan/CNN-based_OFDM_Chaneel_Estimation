function params = configs()
% 一元管理するハイパーパラメータ設定ファイル
params.numTrials = 100; % 100 Reference value for the number of trials when measuring MSE

params.DataSize = 14080; 
params.ValDataSize = 1280; % validation data size

%% CNN structure parameters
params.D = 56;   % 畳み込みチャネル数
params.S = 12;   % サブ機能マップ数
params.M = 4;    % Mapping層数

%% training options
params.MaxEpochs = 10; 
params.MiniBatchSize = 128; % [2^n]の中で最適化
params.InitialLearnRate = 1e-3;
params.L2Regularization = 1e-4;
params.Shuffle ='every-epoch';

params.lossFunction =  "mean-squared-error";



end

%{
% test
params = configs();
disp(params);
%}