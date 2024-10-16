clear; 
tic;

%% variable 
pos=2; % pilot(DM-RS) allocation type of 5G. dmrs-AdditionalPosition= pos1 or pos2
NRB = 20; % subcarrier number = 12*NRB
monte=1;  % 
slots=20; % for ideal & practical LMMSE

%% fixed parameters
scs = 60e3; % subcarrier spacing 
mp = [1:2:12]; % pilot carrier of 1RB
Nfft=512;
% M=4;bps=log2(M);avp=0.5; % bpsk
M=16;bps=log2(M);avp=2.5; % QPSK
Nslot = 1; % number of RBs along time axis
CP = 0.07; %　percentage of cyclic prefic
m_1RB = 12; n_1RB =14; % 1RB= m_1RB carriers * n_1RB time slots
fc=50*10^9; % carrier freqency


matName_FSRCNN=sprintf("../model/FSRCNN_Pos%d_RB%d.mat",pos,NRB);

if pos==1
    np=[3 12]; % pilot timeslot of 1RB
elseif pos==2
    np=[3 8 12];
end

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

%% Monte preparing
SNRdB = 0:5:30;

%% LMMSE prepare
%% also DMRS genaration
dmrsBit=randi([0 1],1,bps*(Ndmrs/length(NP)));
dmrsSym_m=modu(dmrsBit,M);
dmrsSym=[]; 
for i=1:length(NP) dmrsSym=[dmrsSym, dmrsSym_m]; end
% DMRS_LMMSE=diag matrix of 2
dmrsSym_n=repmat(dmrsSym_m(1),length(NP),1);
dmrsDiag_n=diag(dmrsSym_n);
DMRS_LMMSE_n=(dmrsDiag_n*dmrsDiag_n');
dmrsDiag_m=diag(dmrsSym_m);
DMRS_LMMSE_m=(dmrsDiag_m*dmrsDiag_m');




for k=1:length(SNRdB)
    if  SNRdB(k)>=25; Monte=30*monte;
    elseif SNRdB(k)>=15; Monte=20*monte;
    elseif SNRdB(k)>=10 Monte=6*monte;
    else Monte=3*monte;
    end;


    fprintf('%d[dB] %dMonte\n',SNRdB(k),Monte);

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
        maxH=max(max(abs(H_perfect))); 

        normalized_H_perfect=H_perfect/maxH;

        %% LS
        H_LS = rRGseq(dmrs_loc)./dmrsSym; 
        H_LS = reshape(H_LS,length(MP),length(NP));
        
        %% linear
        Hf_linear=zeros(m_1user,length(NP));
        % f axis 
        for i=1:length(NP)
            Hf_linear(:,i) = ...
              interp1(MP,H_LS(:,i),1:m_1user,'linear','extrap');
        end
        % t axis
        H_linear=zeros(m_1user,n_1user);
        for i=1:m_1user
            H_linear(i,:) = ...
              interp1(NP,Hf_linear(i,:),[1:n_1user],'linear','extrap');
        end
        normalized_H_linear=H_linear/maxH;
        H_e = H_perfect-H_linear;
        H_e = H_e(:); H_e(dmrs_loc)=[];
        mse(1,j) = mean(abs(H_e).^2)  ;


        %% perfect LMMSE 
        mode="perfect";
        % get Rh_LS over some slots
        Rh_LS = channel_autocorrelation(mode,t0,DelayProfile,DS,fd,slots,Tofdm,...
        n_1user,Pn,Nofdm,Ng,m_1user,MP,NP, ...
        dmrs_loc,...
        dmrsSym,sOFDM,...
        Nfft);
        % frequency Rh
        I=eye(length(MP));
        H_lmmse_perfect = Rh_LS/(Rh_LS+(Pn*I)/DMRS_LMMSE_m )*H_LS;%
        % QPSK
        % H_lmmse_perfect = Rh_LS/(Rh_LS+2*I*Pn*I)*H_LS;%

        % f axis 
        H_LMMSEf_perfect=zeros(m_1user,length(NP));
        for i=1:length(NP)
            H_LMMSEf_perfect(:,i) = ...
              interp1(MP,H_lmmse_perfect(:,i),1:m_1user,'linear','extrap');
        end
        % t axis
        H_LMMSE_perfect=zeros(m_1user,n_1user);
        for i=1:m_1user
            H_LMMSE_perfect(i,:) = ...
              interp1(NP,H_LMMSEf_perfect(i,:),[1:n_1user],'linear','extrap');
        end

        normalized_H_LMMSE_perfect=H_LMMSE_perfect/maxH;
        H_e = H_perfect-H_LMMSE_perfect;
        H_e = H_e(:); H_e(dmrs_loc)=[];
        mse(2,j) = mean(abs(H_e).^2)  ;

        %% practical LMMSE 
        mode="practical";
        % get Rh_LS over some slots
        Rh_LS = channel_autocorrelation(mode,t0,DelayProfile,DS,fd,slots,Tofdm,...
        n_1user,Pn,Nofdm,Ng,m_1user,MP,NP, ...
        dmrs_loc,...
        dmrsSym,sOFDM,...
        Nfft);
        % frequency Rh
        I=eye(length(MP));
        H_lmmse_practical = Rh_LS/(Rh_LS+(Pn*I)/DMRS_LMMSE_m)*H_LS;%

        % f axis 
        H_LMMSEf_practical=zeros(m_1user,length(NP));
        for i=1:length(NP)
            H_LMMSEf_practical(:,i) = ...
              interp1(MP,H_lmmse_practical(:,i),1:m_1user,'linear','extrap');
        end
        % t axis
        H_LMMSE_practical=zeros(m_1user,n_1user);
        for i=1:m_1user
            H_LMMSE_practical(i,:) = ...
              interp1(NP,H_LMMSEf_practical(i,:),[1:n_1user],'linear','extrap');
        end

        normalized_H_LMMSE_practical=H_LMMSE_practical/maxH;
        H_e = H_perfect-H_LMMSE_practical;
        H_e = H_e(:); H_e(dmrs_loc)=[];
        mse(3,j) = mean(abs(H_e).^2)  ;

 
        %% FSRCNN
        load(matName_FSRCNN);
        nnInput = cat(4,real(H_LS),imag(H_LS));
        H_cnn = predict(channelEstimationCNN,nnInput);
        H_FSRCNN(1:m_1user,1:n_1user) = H_cnn(:,:,1,1)+ 1i*H_cnn(:,:,1,2);
        normalized_H_FSRCNN=H_FSRCNN/maxH;
        H_e = H_perfect-H_FSRCNN;
        H_e=H_e(:); H_e(dmrs_loc)=[];
        mse(4,j) = mean(abs(H_e).^2)  ;

    end
    MSE(:,k)=mean(mse,2);
end     

%% MSE


figure;
markersize=15;
semilogy(SNRdB,MSE(1,:),'r+-','MarkerSize',markersize);hold on;
semilogy(SNRdB,MSE(3,:),'bo-','MarkerSize',markersize);hold on;
semilogy(SNRdB,MSE(2,:),'g^-','MarkerSize',markersize);hold on;grid on;
semilogy(SNRdB,MSE(4,:),'ksquare-','MarkerSize',markersize);hold on;

xlabel('SNR[dB]') 
ylabel('MSE')
legend('LS',...
    'practical LMMSE','ideal LMMSE',...
    '深層学習','FontSize',22);
xlim([0 SNRdB(end)])
set(gca,'FontSize',22)

% channel matrix surface figure
figure;surf(real(H_perfect(:,:)));title('perfect')
figure;surf(real(H_linear(:,:)));title('LS');
figure;surf(real(H_LMMSE_perfect(:,:)));title('ideal LMMSE');
figure;surf(real(H_FSRCNN(:,:)));title('深層学習');


toc;