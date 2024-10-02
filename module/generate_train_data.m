function [trainData,trainLabels,MP] = ...
    generate_train_data(dataSize,np,NRB)
% generate_train_data Generates training data for CNN-based OFDM channel estimation.
%
% Syntax:
%   [trainData, trainLabels, MP] = generate_train_data(dataSize, np, NRB)
%
% Inputs:
%   dataSize - Number of training samples to generate.
%   np       - Pilot positions along the time axis.
%   NRB      - Number of resource blocks along the frequency axis.
%
% Outputs:
%   trainData   - Generated training data.
%   trainLabels - Corresponding labels for the training data.
%   MP          - Pilot positions along the frequency axis.
%
% Example:
%   dataSize = 12800;
%   np = [3, 8, 12];
%   NRB = 20;
%   [trainData, trainLabels, MP] = generate_train_data(dataSize, np, NRB);
%
% Description:
%   This function generates training data for CNN-based OFDM channel estimation.
%   It simulates the OFDM transmission, adds noise, and performs channel estimation
%   to create the training data and labels. The function also calculates the pilot
%   positions along the frequency axis.


%% variable 
scs = 60e3; % 60kHz subcarrier spacing
mp = [1:2:12]; %dmrs position along the frequency axis per RB(resource block)
Nfft=512; % FFT size

%% fixed parameters
M=4;bps=log2(M);avp=0.5; % QPSK
Nslot = 1; % the number of RB along the time axis 
CP = 0.07; % cyclic prefix
m_1RB = 12; n_1RB =14; % RB size=12*14
fc=50*10^9; % carrier frequency

%% Caluculate
Tofdm = (1+CP)/scs; % OFDM symbol duration
m_1user= m_1RB*NRB; % the number of subcarriers
n_1user=n_1RB*Nslot; % the number of OFDM symbols
Ndmrs=NRB*length(mp)* Nslot*length(np); % the number of dmrs symbol
Ns=m_1user*n_1user; % the number of subcarriers*OFDM symbols
Ng=ceil(Nfft*0.07);Nofdm=Nfft+Ng; % guard interval

%% dmrs location
MP=[];NP=[];
for i=1:NRB MP=[MP,(i-1)*m_1RB+mp]; end
for i=1:Nslot NP=[NP,(i-1)*n_1RB+np]; end
a = repmat((NP-1)*m_1user,length(MP),1);
b = repmat(MP',1,length(NP));
c = a+b;
dmrs_loc = c(:);

%% dmrs genaration
dmrsBit=randi([0 1],1,bps*(Ndmrs/length(NP)));
dmrsSym_m=modu(dmrsBit,M);
dmrsSym=[]; 
for i=1:length(NP) dmrsSym=[dmrsSym, dmrsSym_m]; end


trainData=zeros(length(MP),length(NP),2,dataSize);
trainLabels=zeros(m_1user,n_1user,2,dataSize);
for j=1:dataSize
    %% Channel parameters
    SNRdB =  randi([0 30],1,1); 
    t0=randi([0 1e5],1,1); 
    v=randi([0 60],1,1); % km/h
    DelayProfile = char(randsample(["TDL-A","TDL-B","TDL-C",...
      ], 1));
    DS=randsample([55 228],1,1); % delay spread

    fd = v*1000/3600 * fc/(3*10^8); % Doppler frequency

    %% RG
    dataBit = randi([0 1],1,bps*(Ns-Ndmrs));dataSym = modu(dataBit,M);
    [RG,RGseq]=GenRG(Nfft,m_1RB,n_1RB, ...
        NRB,Nslot,...
        dataSym,dmrsSym,dmrs_loc);
    
    %% IFFT,add GI
    IFFTout = ifft(RG);GI=IFFTout((end-Ng+1):end,:);
    sOFDM=([GI; IFFTout]);
    
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
    Pn = 10^(-SNRdB/10)/Nfft *avp;% noise power
    chOut = chOut + awgn*sqrt(Pn/2);
    
    %% remove GI, FFT
    rxSPCOut = reshape(chOut,Nofdm,n_1user);
    rxSPCOut (1:Ng,:) = []; 
    FFTout = fft(rxSPCOut);
    FFTout = FFTout(1:m_1user,:);
    rRGseq = transpose(FFTout(:));
    
    %% perfect channel
    H_perfect=fft(h,Nfft,1); 
    H_perfect=H_perfect(1:m_1user,:);

    %% LS channel estimation
    H_LS = rRGseq(dmrs_loc)./dmrsSym; 
    H_LS = reshape(H_LS,length(MP),length(NP));

    trainData(:,:,1,j)=real(H_LS);
    trainData(:,:,2,j)=imag(H_LS);
    trainLabels(:,:,1,j)=real(H_perfect);
    trainLabels(:,:,2,j)=imag(H_perfect);
end

    
    

