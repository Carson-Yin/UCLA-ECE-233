clc;
clear all;
close all;
N_iter = 1000;%��ͬ������¼������

%�շ������߸���
% Ns = 4;  % �������������
Nt = 128; % �����������
Mt = 4;  % �������Ƶ����
Nr = 32; % ���ն�������
Mr = Mt; % ���ն���Ƶ����

% SNR in dB 
SNR_set = -15:5:15;

% �ŵ�������ɢ��ظ�����ÿ���ذ�������·����
Nc = 5; 
Np = 10; 
sig = deg2rad(10); % ÿ��ɢ����ϽǶȷ���������˹�ֲ��ı�׼��Ƕ�ת��Ϊ���ȣ�

% ��������������Ӧʸ�������ߵ�Ԫ���Ϊ�������
a = @(phi,N) exp(-j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);

BER_Full_Dig = zeros(length(SNR_set),1); 
BER_Hybrid_omp = zeros(length(SNR_set),1);
% BER_Hybrid_rf = zeros(length(SNR_set),1);
% BER_Beem_Steering = zeros(length(SNR_set),1);
for loop=[1]
    Ns=2^loop;
    for SNR_index = 1:length(SNR_set)
        SNR = SNR_set(SNR_index); % in dB 
        SNR_L = 10^(SNR/10); % SNR (���������)

        %Ƶ��Ч�ʳ�ʼ��
%         Temp_ber_full_dig = zeros(length(N_iter),1);
%         Temp_ber_omp = zeros(length(N_iter),1);
%         Temp_ber_rf = zeros(length(N_iter),1); 
    %     Temp_ber_bs = zeros(length(N_iter),1);

        ber_full_dig = 0;
        ber_omp = 0;
%         ber_rf = 0;
%         ber_bs = 0;

        for int=1:N_iter 
            %���ײ��ŵ�
            [H, Phi_AOD, Phi_AOA, Alpha] =...
                            mm_wave_channel(Nt, Nr, Nc, Np, sig);
            atv = Phi_AOD(:);     %������뿪�Ǿ���
            arv = Phi_AOA(:);     %���ն˵���Ǿ���
            A_t = zeros(Nt,Nc*Np);
            A_r = zeros(Nr,Nc*Np);
            for m = 1:Nc*Np
                A_t(:,m) = a(atv(m), Nt);%�����������Ӧʸ������
                A_r(:,m) = a(arv(m), Nr);%���ն�������Ӧʸ������
            end        
            %��С���������ջ�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            % ����ȫ����Ԥ���롢�����
            [U,S,V]=svd(H); 
            F_opt = V(:,1:Ns);
            W_opt = U(:,1:Ns);
    %         W_opt = ((1/sqrt(SNR_L))*(F_opt'*H'*H*F_opt+Ns/SNR_L*eye(Ns))\(F_opt'*H'))';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %Beem Steeringģ��Ԥ����
%             [F_BS,W_BS] = Beamsteering(H,A_t,A_r,Ns);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % OMP���Ԥ���롢�����
            [F_R_omp, F_B_omp] = omp(F_opt, A_t, Ns, Mt);
            CovRx_omp = (SNR_L/Ns)*H*F_R_omp*F_B_omp*F_B_omp'*F_R_omp'*H'+eye(Nr);
            W_MMSE_omp = ((1/sqrt(SNR_L))*(F_B_omp'*F_R_omp'*H'*H*F_R_omp*F_B_omp+Ns/SNR_L*eye(Ns))\(F_B_omp'*F_R_omp'*H'))';
            [W_R_omp, W_B_omp] = MMSE_omp(W_MMSE_omp, A_r, CovRx_omp, Mr);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % RFԤ���롢�����
%             [F_R_rf, F_B_rf] = RF(F_opt, Ns, Mt, Nt); 
%             CovRx_rf = (SNR_L/Ns)*H*F_R_rf*F_B_rf*F_B_rf'*F_R_rf'*H'+eye(Nr);
%             W_MMSE_rf = ((1/sqrt(SNR_L))*(F_B_rf'*F_R_rf'*H'*H*F_R_rf*F_B_rf+Ns/SNR_L*eye(Ns))\(F_B_rf'*F_R_rf'*H'))';
%             [W_R_rf, W_B_rf] = MMSE_RF(W_MMSE_rf, Ns, Mr, Nr, CovRx_rf); 
            %%%%%%%%%%%%%%%%-----�����-----%%%%%%%%%%%%%%
            %��������
            %%%%%%%%%%%%%%%%-----BPSK-----%%%%%%%%%%%%%%      
    %         S1=rand(Ns,1);
    %         S=round(S1);
    %         data=2.*S-1;
            %%%%%%%%%%%%%%%%-----4QAM-----%%%%%%%%%%%%%%
            msg = randsrc(Ns, 1, (0:3));
            data = qammod(msg, 4)/(sqrt(2));
            %Ԥ����
            X_opt=F_opt*data;%��վ�˵ķ����ź�
            X_omp=F_R_omp*F_B_omp*data;
%             X_rf=F_R_rf*F_B_rf*data;
%             X_bs=F_BS*data;
            %ͨ���ŵ�����
            Y_opt = awgn(sqrt(1/Ns)*H*X_opt, SNR, 'measured');
            Y_omp = awgn(sqrt(1/Ns)*H*X_omp, SNR, 'measured');
%             Y_rf = sqrt(P/Ns)*H*X_rf+Noise;
%             Y_bs = sqrt(P)*H*X_bs+Noise;
            S_opt = W_opt'*Y_opt;
            S_omp = W_B_omp'*W_R_omp'*Y_omp;
%             S_rf = W_B_rf'*W_R_rf'*Y_rf;
%             S_bs = W_BS'*Y_bs;
            %���ն�����
    %         [Temp_ber_full_dig(i,1)] = ber_bpsk(S_opt, Ns, data);  
    %         [Temp_ber_omp(i,1)] = ber_bpsk(S_omp, Ns, data); 
    %         [Temp_ber_rf(i,1)] = ber_bpsk(S_rf, Ns, data);
    %         [Temp_ber_bs(i,1)] = ber(S_bs, Ns, data);

            [Temp_ber_full_dig(int,1)] = biterr(msg,qamdemod(S_opt, 4));
            [Temp_ber_omp(int,1)] = biterr(msg,qamdemod(S_omp, 4));
%             [Temp_ber_rf(i,1)] = ber_4QAM(S_rf, Ns, data1);
%             [Temp_ber_bs(i,1)] = ber_4QAM(S_bs, Ns, data1);

            ber_full_dig=ber_full_dig+Temp_ber_full_dig(int,1);
            ber_omp=ber_omp+Temp_ber_omp(int,1);
%             ber_rf=ber_rf+Temp_ber_rf(i,1);
%             ber_bs=ber_bs+Temp_ber_bs(i,1);
        end
        BER_Full_Dig(SNR_index) = real(ber_full_dig)/(N_iter*Ns*2);
        BER_Hybrid_omp(SNR_index) = real(ber_omp)/(N_iter*Ns*2);       
%         BER_Hybrid_rf(SNR_index) = real(ber_rf)/(N_iter*Ns*2); 
%         BER_Beem_Steering(SNR_index) = real(ber_bs)/(N_iter*Ns*2);
    end
    figure(2);
    semilogy(SNR_set,BER_Full_Dig,'-', 'Linewidth', 1.5,'MarkerSize',4);
    hold on;
    semilogy(SNR_set,BER_Hybrid_omp,'.-', 'Linewidth', 1.5,'MarkerSize',4);
%     semilogy(SNR_set,BER_Hybrid_rf,'g^:', 'Linewidth', 1.5,'MarkerSize',4);
%     semilogy(SNR_set,BER_Beem_Steering,'bx:', 'Linewidth', 1.5,'MarkerSize',4);
    set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
    set(get(gca,'YLabel'),'String','BER','Interpreter','latex');
    hl = legend('Digital','OMP','Location','Northwest');
    set(hl, 'Fontsize', 12,'Interpreter','latex');
    grid on;
end






