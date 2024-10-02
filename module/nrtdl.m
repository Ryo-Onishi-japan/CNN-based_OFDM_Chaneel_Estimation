% revisited Jakes model
% 4096FFT
function [T_WH] = nrtdl(fd,Tofdm,Ns,t0,DelayProfile,DS)
% Tofdm : sampling time[s]
% Ns : ofdm symbols
% fd : max doppler[Hz]



%% 4096FFT
% DS=55;
switch DS
    case 55
        switch DelayProfile
           case 'TDL-A'
                P_index=[0	1	3	4	5	7	8	9	16]+1;
                PDP=[0.045708819	0.398107171	0.218776162	0.030199517	0.074131024	0.05370318	0.012882496	0.01023293	0.001071519];
                P_tau=zeros(1,P_index(end)); 
                P_tau(P_index)=PDP;
            case 'TDL-B'
                P_index=[0	1	2	3	5	6	7	8]+1;
                PDP=[1	0.128824955	0.331131121	0.645654229	0.104712855	0.072443596	0.032359366	0.074131024];
                P_tau=zeros(1,P_index(end)); 
                P_tau(P_index)=PDP;      
            case 'TDL-C' 
                P_index=[0	1	2	4	5	7	8	9	11	12	15]+1;
                PDP=[0.363078055	1	0.309029543	0.134896288	0.047863009	0.040738028	0.040738028	0.02630268	0.026915348	0.00691831	0.005248075];
                P_tau=zeros(1,P_index(end)); 
                P_tau(P_index)=PDP;      
            % case 'TDL-D'
            %     P_index=[0	1	2	3	4	7	13	16	21]+1;
            %     PDP=[0.999660945	0.007943282	0.005248075	0.005128614	0.006456542	0.001659587	0.004365158	0.003311311	0.001698244];
            %     P_tau=zeros(1,P_index(end)); 
            %     P_tau(P_index)=PDP;      
            % case 'TDL-E'  
            %     P_index=[0	1	3	4	6	9	20	35]+1;
            %     PDP=[0.999382187	0.010471285	0.013803843	0.005888437	0.002754229	0.009549926	0.001047129	0.001202264];
            %     P_tau=zeros(1,P_index(end)); 
            %     P_tau(P_index)=PDP;      
            otherwise
                warning('Unexpected TDL model')
         end 
    case 228
        switch DelayProfile
            case 'TDL-A'
                P_index=[0	3	4	5	11	13	15	16	17	18	21	29	31	32	34	35	37	68
]+1;
                PDP=[0.045708819	0.602559586	0.089125094	0.102329299	0.025703958	0.218776162	0.057543994	0.021379621	0.030199517	0.083176377	0.074131024	0.05370318	0.023988329	0.014791084	0.012882496	0.021877616	0.01023293	0.001071519
];
                P_tau=zeros(1,P_index(end)); 
                P_tau(P_index)=PDP;      
            case 'TDL-B'
                P_index=[0	1	2	3	4	8	9	11	12	14	20	21	25	29	30	34
]+1;
                PDP=[1	0.602559586	0.104712855	0.45708819	0.128824955	0.331131121	0.26915348	0.177827941	0.645654229	0.173780083	0.060255959	0.104712855	0.072443596	0.032359366	0.120226443	0.074131024
];
                P_tau=zeros(1,P_index(end)); 
                P_tau(P_index)=PDP;         
            case 'TDL-C' 
                P_index=[0	1	2	4	5	6	7	9	15	19	30	32	38	39	44	46	49	61
]+1;
                PDP=[0.363078055	0.758577575	0.301995172	1	0.181970086	0.085113804	0.077624712	0.208929613	0.134896288	0.047863009	0.040738028	0.040738028	0.02630268	0.019498446	0.025118864	0.026915348	0.00691831	0.005248075
];
                P_tau=zeros(1,P_index(end)); 
                P_tau(P_index)=PDP;      
%             case 'TDL-D'
%                 P_index=[0	4	10	12	13	18	28	56	66	68	88
% ]+1;
%                 PDP=[0.999660945	0.007943282	0.016218101	0.005128614	0.009772372	0.006456542	0.001659587	0.004365158	0.003311311	0.001	0.001698244
% ];
%                 P_tau=zeros(1,P_index(end)); 
%                 P_tau(P_index)=PDP;      
%             case 'TDL-E'  
%                 P_index=[0	4	5	13	14	19	26	38	84	145
% ]+1;
%                 PDP=[0.999382187	0.010471285	0.005754399	0.013803843	0.005495409	0.005888437	0.002754229	0.009549926	0.001047129	0.001202264
% ];
%                 P_tau=zeros(1,P_index(end)); 
%                 P_tau(P_index)=PDP;      
             otherwise
                warning('Unexpected TDL model')
         end

    otherwise
        warning('Unexpected DS value')
end

% t0 = 0; % Removing randomness
t=t0*ones(1,Ns)+Tofdm*[0:Ns-1];

wd = 2*pi*fd;
Ndelay = length(P_tau);
N0= 2^(ceil( log(Ndelay)/log(2)) );
N=4*N0;   % oscillator wave
WHMatrix = hadamard(N0);
alpha = 2*pi*([1:N0]-0.5)/N; 
beta = pi*[1:N0]/N0;
theta = rand(N0,1)*2*pi;
% theta = zeros(N0,1); % Removing randomness
w=wd*cos(alpha);

complex = cos(beta)+1i*sin(beta);
in_cos = transpose(w)*t+theta*ones(1,Ns);
in_sigma = transpose(complex).*cos(in_cos);
%T = sqrt(2/N0)*sum(in_sigma); %single path。あとRayleighの期待値は0.886

T_WH = zeros(Ndelay,Ns);
for i=1:Ndelay
    in_sigma_WH = WHMatrix(i,:)'.*in_sigma;
    %in_sigma_WH = WHMatrix(:,i).*in_sigma;
    T_WH(i,:) = sqrt(P_tau(i))*sqrt(2/N0)*sum(in_sigma_WH);
end

% clear,close
% fd = 926;
% Ts = 1e-6;
% Ns = 50000;
% h = Modified_Jakes(fd,Ts,Ns);
% %figure;plot([1:Ns],10*log10(abs(h)))
% mean(abs(h),2)




