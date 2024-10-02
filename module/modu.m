function x=modu(u,M)
%入出は1行ベクトル   
if M == 2
    %行ベクトルがシンボル
    x = 2*u -1; %0,1→-1,+1
else

    row_count = size(u,1);
    x = [];
    for i = 1:row_count
        u0 = u(i,:);
        bits=log2(M); %1symbolのbit数
        %symbol数
        num_symbols = numel(u0) / bits;
        %行ベクトルがシンボル
        u_bits = reshape(u0, bits, num_symbols)'; % [1 2
        % 　　　　　　　　　　　　　　　　　　　　　　　3 4
        % 　　　　　　　　　　　　　　　　　　　　　　　5 6 
    
        %I,Q相で分離
        u_Re = u_bits(:, 1:(bits/2));
        u_Im = u_bits(:,(bits/2+1):end);
        %uをグレイコードとして扱い、binaryに変換し10進数に変換
    
        % Gray→Binary
        s_ReB(:,1) = u_Re(:,1);
        s_ImB(:,1) = u_Im(:,1);
        for j = 2:bits/2
            s_ReB(:,j) = xor(u_Re(:,j), s_ReB(:,j-1));
            s_ImB(:,j) = xor(u_Im(:,j), s_ImB(:,j-1));
        end
        
        %Binary→Decimal
        s_ReD = bi2de(s_ReB, 'left-msb') ;
        s_ImD = bi2de(s_ImB, 'left-msb') ;
        %信号の中心を原点に 
        s_ReD = s_ReD - (bits/2-0.5);
        s_ImD = s_ImD - (bits/2-0.5);
        %送信信号
        x = [x,(s_ReD + 1i*s_ImD)];
    end  
    %列ベクトル→行ベクトル
    x =transpose(x);
end

end