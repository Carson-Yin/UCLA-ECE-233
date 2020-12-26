function [H_n, Phi_AOD, Phi_AOA, Alpha] = mm_wave_channel(Nt, Nr, Nc, Np, sig)
a = @(phi,N) exp(-j*pi*sin(phi)*(0:1:N-1)).'/sqrt(N);%均匀线性阵列
% The azimuth angles for the cluster's centers
Phi_AOD_m = 2*pi*rand(Nc,1);                  %离开角均值（对于每一个散射簇）
Phi_AOA_m = 2*pi*rand(Nc,1);                  %到达角均值（对于每一个散射簇）
Alpha = (1/sqrt(2))*(randn(Nc,Np)+j*randn(Nc,Np));%信道增益因子
% 初始化
Phi_AOD = zeros(Nc,Np);
Phi_AOA = zeros(Nc,Np);
H = zeros(Nr,Nt);
 for i=1:Nc
     % 每个散射簇不同散射路径的离开角和到达角    
     phi_AOD = tr_laprnd(Np, 1, Phi_AOD_m(i), sig);
     phi_AOA = tr_laprnd(Np, 1, Phi_AOA_m(i), sig);    
         for l=1:Np
         Phi_AOD(i,l) =  phi_AOD(l);
         Phi_AOA(i,l) =  phi_AOA(l);
         alpha = Alpha(i,l);
         A_t = a(phi_AOD(l), Nt); % 发射端天线阵列响应矢量
         A_r = a(phi_AOA(l), Nr); % 接收端天线阵列响应矢量
         H = H + alpha*A_r*A_t'; % 信道矩阵
         end
 end
 H_n = sqrt(Nt*Nr/(Nc*Np))*H;    %归一化信道矩阵
end

function y = tr_laprnd(m, n, mu, sigma)
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b*sign(u).*log(1- 2*(1-exp(-pi/b))*abs(u));
end