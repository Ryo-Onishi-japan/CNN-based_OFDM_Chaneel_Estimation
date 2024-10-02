function [Rh_LS] = ...
    Tracking_H_LS(mode,t0,DelayProfile,DS,fd,slots,Tofdm,...
    n_1user,Pn,Nofdm,Ng,m_1user,MP,NP,...
    dmrs_loc,...
    dmrsSym,sOFDM,...
    Nfft)

%% fixed
n_1RB =14;

%% caluculation
duration_slot=n_1RB*Tofdm;
t=t0-slots*Tofdm;

for j=1:slots
    t=t+(j-1)*Tofdm;   
    
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
    h = H_perfect(:);
    hp=h(dmrs_loc); 
    hp=reshape(hp,length(MP),length(NP));

    if mode =='perfect'
        % RHfp
        for i=1:length(NP)
            index=i+(j-1)*size(hp,2);
            rhfp(:,:,index) = hp(:, i)* hp(:, i)';
        end
    elseif mode =='practical'
        % LS
        H_LS = rRGseq(dmrs_loc)./dmrsSym; 
        H_LS = reshape(H_LS,length(MP),length(NP));
        for i=1:size(H_LS,2);
            index=i+(j-1)*size(H_LS,2);
            rhfp(:,:,index)= H_LS(:, i)* H_LS(:, i)';
        end
    end
end

Rh_LS=mean(rhfp,3);