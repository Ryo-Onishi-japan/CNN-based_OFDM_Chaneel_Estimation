function [trainData,trainLabels,MP] = ...
    TrainData_FSRCNN(dataSize,np,NRB)

%% variable 
scs = 60e3;
mp = [1:2:12];
Nfft=512;

%% fixed parameters
M=4;bps=log2(M);avp=0.5;
Nslot = 1; % RG
CP = 0.07; % [%]
m_1RB = 12; n_1RB =14;
fc=50*10^9; 


%% Caluculate
Tofdm = (1+CP)/scs; % OFDM symbol duration
m_1user= m_1RB*NRB; n_1user=n_1RB*Nslot;
Ndmrs=NRB*length(mp)* Nslot*length(np); 
Ns=m_1user*n_1user; % data+dmrs symbol
% Number of samples transmitted by base station
Ng=ceil(Nfft*0.07);Nofdm=Nfft+Ng; % symbols/ofdm symbol

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
    %% randomness
    SNRdB =  randi([0 30],1,1); 
    t0=randi([0 1e5],1,1); 
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
    Pn = 10^(-SNRdB/10)/Nfft *avp;%/Rofdm
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

    trainData(:,:,1,j)=real(H_LS);
    trainData(:,:,2,j)=imag(H_LS);
    trainLabels(:,:,1,j)=real(H_perfect);
    trainLabels(:,:,2,j)=imag(H_perfect);
end

    
    

