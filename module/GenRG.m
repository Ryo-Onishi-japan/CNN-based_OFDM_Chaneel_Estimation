% (frequency,time)=(m,n)
function[RG,RGseq]=GenRG(Nfft, ...
    m_1RB,n_1RB, ... % 1RB size definition
    NRB,Nslot,... 
    dataSym,dmrsSym,dmrs_loc)

m_1user=m_1RB*NRB;n_1user=n_1RB*Nslot;
Ns = length(dataSym)+length(dmrsSym);
%% RG
RGseq_1user=zeros(1,Ns);
idata=1;idmrs=1;
for i = 1:Ns
    if i>dmrs_loc(end)
        
        RGseq_1user(i)=dataSym(idata);
        idata=idata+1;
    elseif i == dmrs_loc(idmrs)
        RGseq_1user(i)=dmrsSym(idmrs);
        idmrs = idmrs+1;
        
    else
        RGseq_1user(i)=dataSym(idata);
        idata=idata+1;
    end
end
RG_1user = reshape(RGseq_1user,m_1user,n_1user);
RG=zeros(Nfft,n_1user);
RG(1:m_1user,:)=RG_1user;
RGseq=reshape(RG,1,Nfft*Nslot*n_1RB);

% figure;
% surf(abs(RG)); 
% xlabel('OFDM symbol','FontSize',20)%120404
% ylabel('Frequency','FontSize',20)%120404
% title('RBs of 1user')
% set(gca,'FontSize',18)

end

