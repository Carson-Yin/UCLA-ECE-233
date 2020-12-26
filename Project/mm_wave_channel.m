function [H_n, Phi_AOD, Phi_AOA, Alpha] = mm_wave_channel(Nt, Nr, Nc, Np, sig)
a = @(phi,N) exp(-j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);%������������
% The azimuth angles for the cluster's centers
Phi_AOD_m = 2*pi*rand(Nc,1);                  %�뿪�Ǿ�ֵ������ÿһ��ɢ��أ�
Phi_AOA_m = 2*pi*rand(Nc,1);                  %����Ǿ�ֵ������ÿһ��ɢ��أ�
Alpha = (1/sqrt(2))*(randn(Nc,Np)+j*randn(Nc,Np));%�ŵ���������
% ��ʼ��
Phi_AOD = zeros(Nc,Np);
Phi_AOA = zeros(Nc,Np);
H = zeros(Nr,Nt);
 for i=1:Nc
     % ÿ��ɢ��ز�ͬɢ��·�����뿪�Ǻ͵����    
     phi_AOD = tr_laprnd(Np, 1, Phi_AOD_m(i), sig);
     phi_AOA = tr_laprnd(Np, 1, Phi_AOA_m(i), sig);    
         for l=1:Np
         Phi_AOD(i,l) =  phi_AOD(l);
         Phi_AOA(i,l) =  phi_AOA(l);
         alpha = Alpha(i,l);
         A_t = a(phi_AOD(l), Nt); % ���������������Ӧʸ��
         A_r = a(phi_AOA(l), Nr); % ���ն�����������Ӧʸ��
         H = H + alpha*A_r*A_t'; % �ŵ�����
         end
 end
 H_n = sqrt(Nt*Nr/(Nc*Np))*H;    %��һ���ŵ�����
end

function y = tr_laprnd(m, n, mu, sigma)
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b*sign(u).*log(1- 2*(1-exp(-pi/b))*abs(u));
end