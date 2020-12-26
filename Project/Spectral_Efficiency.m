clc;
clear all;
close all;
N_iter = 1000;%相同信噪比下计算次数

%收发端天线个数
Ns = 2;  % 发射端数据流数
Nt = 128; % 发射端天线数
Mt = 8;  % 发射端射频链数
Nr = 32; % 接收端天线数
Mr = Mt; % 接收端射频链数

% SNR in dB 
SNR_set = -30:0.1:5;

% 信道特征：散射簇个数和每个簇包含传播路径数
Nc = 5; 
Np = 10; 
sig = deg2rad(10); % 每个散射簇上角度符合拉普拉斯分布的标准差（角度转换为弧度）

% 均匀线性阵列响应矢量，天线单元间距为半个波长
a = @(phi,N) exp(-j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);

% 频谱效率初始化
SE_Full_Dig = zeros(length(SNR_set),1); 
SE_Hybrid_omp = zeros(length(SNR_set),1);
% SE_Hybrid_rf = zeros(length(SNR_set),1);
% SE_Beem_Steering = zeros(length(SNR_set),1);
    
for SNR_index = 1:length(SNR_set)
    SNR = SNR_set(SNR_index); % in dB 
    SNR_L = 10^(SNR/10); % SNR (线性信噪比)
    
    %频谱效率初始化
    Temp_se_full_dig = 0;
    Temp_se_omp = 0;
%     Temp_se_rf = 0;
%     Temp_se_bs = 0;
    
    for i=1:N_iter 
        %毫米波信道
        [H, Phi_AOD, Phi_AOA, Alpha] =...
                        mm_wave_channel(Nt, Nr, Nc, Np, sig);
        atv = Phi_AOD(:);     %发射端离开角矩阵
        arv = Phi_AOA(:);     %接收端到达角矩阵
        A_t = zeros(Nt,Nc*Np);
        A_r = zeros(Nr,Nc*Np);
        for m = 1:Nc*Np
            A_t(:,m) = a(atv(m), Nt);%发射端阵列响应矢量矩阵
            A_r(:,m) = a(arv(m), Nr);%接收端阵列响应矢量矩阵
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % 最优全数字预编码
        [U,S,V]=svd(H); 
        F_opt = V(:,1:Ns);
        W_opt = ((1/sqrt(SNR_L))*(F_opt'*H'*H*F_opt+Ns/SNR_L*eye(Ns))\(F_opt'*H'))';
        Rn_opt = W_opt'*W_opt;
        R_dig = log2(det(eye(Ns)+(SNR_L/Ns)*(Rn_opt\(W_opt'*H*F_opt*F_opt'*H'*W_opt))));
        %互信息量（频谱效率）
        Temp_se_full_dig = Temp_se_full_dig + R_dig;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Beem Steering模拟预编码
%         [F_BS,W_BS] = findSteeringVector(H,A_t,A_r,Ns);
%         Rn_bs = W_BS'*W_BS;
%         R_bs = log2(det(eye(Ns)+(SNR_L/Ns)*(Rn_bs\(W_BS'*H*F_BS*F_BS'*H'*W_BS))));
%         Temp_se_bs = Temp_se_bs + R_bs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OMP混合预编码
        [F_R_omp, F_B_omp] = omp(F_opt, A_t, Ns, Mt);
        % OMP混合组合器
        CovRx_omp = (SNR_L/Ns)*H*F_R_omp*F_B_omp*F_B_omp'*F_R_omp'*H'+eye(Nr);
        W_MMSE_omp = ((1/sqrt(SNR_L))*(F_B_omp'*F_R_omp'*H'*H*F_R_omp*F_B_omp+Ns/SNR_L*eye(Ns))\(F_B_omp'*F_R_omp'*H'))';
        [W_R_omp, W_B_omp] = MMSE_omp(W_MMSE_omp, A_r, CovRx_omp, Mr);
        Rn_omp = W_B_omp'*W_R_omp'*W_R_omp*W_B_omp;
        R_omp = log2(det(eye(Ns)+(SNR_L/Ns)*(Rn_omp\(W_B_omp'*W_R_omp'*H*F_R_omp*F_B_omp*F_B_omp'*F_R_omp'*H'*W_R_omp*W_B_omp))));
        Temp_se_omp = Temp_se_omp + R_omp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%         % RF预编码
%         [F_R_rf, F_B_rf] = RF(F_opt, Ns, Mt, Nt); 
%         % RF混合组合器
%         CovRx_rf = (SNR_L/Ns)*H*F_R_rf*F_B_rf*F_B_rf'*F_R_rf'*H'+eye(Nr);
%         W_MMSE_rf = ((1/sqrt(SNR_L))*(F_B_rf'*F_R_rf'*H'*H*F_R_rf*F_B_rf+Ns/SNR_L*eye(Ns))\(F_B_rf'*F_R_rf'*H'))';
%         [W_R_rf, W_B_rf] = MMSE_RF(W_MMSE_rf, Ns, Mr, Nr, CovRx_rf); 
%         Rn_rf = W_B_rf'*W_R_rf'*W_R_rf*W_B_rf;
%         R_rf = log2(det(eye(Ns)+(SNR_L/Ns)*(Rn_rf\(W_B_rf'*W_R_rf'*H*F_R_rf*F_B_rf*F_B_rf'*F_R_rf'*H'*W_R_rf*W_B_rf))));
%         Temp_se_rf = Temp_se_rf + R_rf;  
    end  
    SE_Full_Dig(SNR_index) = real(Temp_se_full_dig)/N_iter;
    SE_Hybrid_omp(SNR_index) = real(Temp_se_omp)/N_iter;       
%     SE_Hybrid_rf(SNR_index) = real(Temp_se_rf)/N_iter; 
%     SE_Beem_Steering(SNR_index) = real(Temp_se_bs)/N_iter; 
end
figure(1);
plot(SNR_set,SE_Full_Dig,'-', 'Linewidth', 1.5,'MarkerSize',4);
hold on;
plot(SNR_set,SE_Hybrid_omp,'.-', 'Linewidth', 1.5,'MarkerSize',4);
% plot(SNR_set,SE_Hybrid_rf,'g^:', 'Linewidth', 1.5,'MarkerSize',4);
% plot(SNR_set,SE_Beem_Steering,'bx:', 'Linewidth', 1.5,'MarkerSize',4);
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Spectral Efficiency (bps/Hz)','Interpreter','latex');
hl = legend('Optimal','OMP','Location','Northwest');
set(hl, 'Fontsize', 12,'Interpreter','latex');
grid on;


