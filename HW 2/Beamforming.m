clc;
clear all;
% close all;
% N_T = 16;        %Number of transmit antennas
% f_c = 2.4e9;              %Carrier frequency (Hz)
% lambda = 3e8/f_c;         %Wavelength
% d = [lambda/2, lambda*2, lambda/8];             %Antenna spacing
% phi = -pi/6;                  %Steering angle (radians)
% theta = -pi/2:0.01:pi/2;  %Angles that we will use to compute beamforming gains
% figure;
% for n=1:length(d)
%     gain = zeros(1,length(theta));                                                      %Initialize gain vector   
%     w_BF = exp(-1j*2*pi*(0:N_T-1)'*d(n)*sin(phi)/lambda);                               %Beamforming vector
%     for t = 1:length(theta)
%         a_theta = exp(-1j*2*pi*(0:N_T-1)'*d(n)*sin(theta(t))/lambda)/sqrt(N_T);      %ULA array response  
%         gain(t) = abs(w_BF'*a_theta)^2;                                                 %Array gain
%     end
%     gain = gain./(max(gain));
%     gain = 10*log10(gain);
%     polarplot(theta, gain)
%     hold on
% end
% rlim('auto')
% thetalim([-90 90])
% rlim([-10 0])
% legend('\lambda/2','2\lambda','\lambda/8')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nt = 32;
% L = 4;
% Pt = 1;
% B = 200000;
% N0 = 4e-21;
% Cr = zeros(1,15);
% Cs = zeros(1,15);
% Cr2 = zeros(1,15);
% Cs2 = zeros(1,15);
% Cr3 = zeros(1,15);
% Cs3 = zeros(1,15);
% Ftr = zeros(Nt,Nt);
% for i = 1:Nt
%     Ftr(:,i) = sqrt(1/Nt)* tr_laprnd(-pi/2+(i-1)*pi/Nt,Nt);
% end
% 
% P = 2;
% for Nr = 4:2:32
%     K = Nr;
%     Wtr = zeros(Nr,Nr);
%     for j = 1:Nr
%         Wtr(:,j) = sqrt(1/Nr)* tr_laprnd(-pi/2+(j-1)*pi/Nr,Nr);
%     end
%     for loop = 1:1000
%         Alpha = normrnd(0,max(Nt,Nr)/2,[1 L])+j*normrnd(0,max(Nt,Nr)/2,[1 L]);
%         Theta = pi*(rand(1,L)-0.5);                  %离开角均值（对于每一个散射簇）
%         Phi = pi*(rand(1,L)-0.5);                  %到达角均值（对于每一个散射簇）
%         Hr = normrnd(0,1/2,[Nr Nt])+j*normrnd(0,1/2,[Nr Nt]);
%         Hs = zeros(Nr,Nt);    
%         for i=1:L 
%             A_t = sqrt(1/Nt)*tr_laprnd(Phi(i), Nt); % 发射端天线阵列响应矢量
%             A_r = sqrt(1/Nr)*tr_laprnd(Theta(i), Nr); % 接收端天线阵列响应矢量
%             Hs = Hs + Alpha(i)*A_r*A_t'; % 信道矩阵
%         end
%         [U1,S1,V1]=svd(Hr);
%         [U2,S2,V2]=svd(Hs);
%         for k = 1:K
%             F1 = V1(:,k);
%             W1 = U1(:,k)';
%             Cr(Nr/2-1) = Cr(Nr/2-1)+log2(1+(Nt*Nr*Pt*(abs(W1*Hr*F1)^2)/(K*B*N0)));
%             F2 = V2(:,k);
%             W2 = U2(:,k)';
%             Cs(Nr/2-1) = Cs(Nr/2-1)+log2(1+(Nt*Nr*Pt*(abs(W2*Hs*F2)^2)/(K*B*N0)));
%         end
%         for ii = 1:Nt
%             for jj = 1:Nr
%                 y1(ii,jj) = abs(Wtr(:,jj)'*Hr*Ftr(:,ii))^2;
%                 y2(ii,jj) = abs(Wtr(:,jj)'*Hs*Ftr(:,ii))^2;
%             end
%         end
%         Y1 = sort(y1(:));
%         Y2 = sort(y2(:));
%         [row11,col11]=find(y1==Y1(end));
%         [row12,col12]=find(y1==Y1(end-1));
%         [row21,col21]=find(y2==Y2(end));
%         [row22,col22]=find(y2==Y2(end-1));
%         U3 = [Wtr(:,col11) Wtr(:,col12)];
%         V3 = [Ftr(:,row11) Ftr(:,row12)];
%         U4 = [Wtr(:,col21) Wtr(:,col22)];
%         V4 = [Ftr(:,row21) Ftr(:,row22)];
%         for p = 1:P
%             F5 = V1(:,p);
%             W5 = U1(:,p)';
%             Cr2(Nr/2-1) = Cr2(Nr/2-1)+log2(1+(Nt*Nr*Pt*(abs(W5*Hr*F5)^2)/(P*B*N0)));
%             F6 = V2(:,p);
%             W6 = U2(:,p)';
%             Cs2(Nr/2-1) = Cs2(Nr/2-1)+log2(1+(Nt*Nr*Pt*(abs(W6*Hs*F6)^2)/(P*B*N0)));
%             
%             F3 = V3(:,p);
%             W3 = U3(:,p)';
%             Cr3(Nr/2-1) = Cr3(Nr/2-1)+log2(1+(Nt*Nr*Pt*(abs(W3*Hr*F3)^2)/(P*B*N0)));
%             F4 = V4(:,p);
%             W4 = U4(:,p)';
%             Cs3(Nr/2-1) = Cs3(Nr/2-1)+log2(1+(Nt*Nr*Pt*(abs(W4*Hs*F4)^2)/(P*B*N0)));
%         end
%     end
% end
% Cr = B*Cr/loop;
% Cs = B*Cs/loop;
% Cr2 = B*Cr2/loop;
% Cs2 = B*Cs2/loop;
% Cr3 = B*Cr3/loop;
% Cs3 = B*Cs3/loop;
% 
% figure
% Nr = 4:2:32;
% semilogy(Nr,Cr,'-bo',Nr,Cs,'-b*');
% % plot(Nr,Cr,Nr,Cs);
% title('Channels capacities v.s. number of receive antennas');
% xlabel('N_r');
% ylabel('C_r');
% legend('H_r channrl capacity','H_s channrl capacity')
% 
% figure
% Nr = 4:2:32;
% semilogy(Nr,Cr,'-bo',Nr,Cs,'-b*',Nr,Cr2,'-ro',Nr,Cs2,'-r*',Nr,Cr3,'-go',Nr,Cs3,'-g*');
% % plot(Nr,Cr,Nr,Cs);
% title('Channels capacities and achievable data rate');
% xlabel('N_r');
% ylabel('C_r');
% legend('H_r channrl capacity','H_s channrl capacity','H_r achievable rate (known)','H_s achievable rate (known)','H_r achievable rate (unknown)','H_r achievable rate (unknown)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angle=-90:1:90;
theta = angle*pi/180;
g = zeros(4,length(theta));
psd = zeros(4,4096);
count = 1;
b = 1;
% beta_3 = 0;
for M = [8 32];
%     for b = [12 4];
    for beta_3 = [0 -133];
        clear S;
        clear X;
        clear P;
        clear H;
        beta_1 = 1;

        
        s1 = qammod(randsrc(1, 1000, (0:3)), 4)/(sqrt(2));
%         s2 = qammod(randsrc(1, 1000, (0:3)), 4)/(sqrt(2));
        s2 = zeros(1,1000);
        transmit_filter = rcosdesign(0.5,8,4,'sqrt');
        s1_t = upsample(s1,4);
        s2_t = upsample(s2,4);
        s1_tilde = conv(s1_t,transmit_filter);
        s2_tilde = conv(s2_t,transmit_filter);
        h1 = tr_laprnd(pi/6,M);
        h2 = tr_laprnd(2*pi/9,M);
        H = [h1,h2]';
        S = [s1_tilde;s2_tilde];
        P = H'*(inv(H*(H')));
        X = P*S;
        delta = (2/M)/2^b;
        codebook = -1/M-delta/2:delta:1/M+delta/2;

        for antenna = 1:M
            [~,Y1(antenna,:)] = quantiz(real(X(antenna,:)),linspace(-1/M, 1/M, 2^b),linspace(-1/M, 1/M, 2^b+1));
            [~,Y2(antenna,:)] = quantiz(imag(X(antenna,:)),linspace(-1/M, 1/M, 2^b),linspace(-1/M, 1/M, 2^b+1));
            
        end
%         Y = Y1 + j* Y2;
        Y = X;
        z = beta_1.*Y+beta_3.*Y.*abs(Y).^2;

        [Pxx,w]=periodogram(z(1,:),'PSD');
        psd(count,:) = Pxx.';
   
        for ang = 1:length(theta)
            a = tr_laprnd(theta(ang),M)';
            for i = 1:4032 
                g(count,ang) = g(count,ang)+ (1/4032).*(abs((a*z(:,i))^2));
            end
        end
        count = count+1;
    end
    
end
figure
for count =1:4
    psd(count,:) = 10*log10(psd(count,:));
    plot(w/pi,psd(count,:));
    xlabel('w/\pi rad');
    ylabel('PSD/dB');
    title('Power spectral density of the signal');
    axis([0 2 -150 0])
    grid on
    hold on
end
% legend('M=8,b=12','M=8,b=4','M=32,b=12','M=32,b=4')
legend('M=8,\beta_3=0','M=8,\beta_3=-133','M=32,\beta_3=0','M=32,\beta_3=-133')
figure
for count =1:4
    gain(count,:) = 10*log10(g(count,:));
    plot(theta/pi*180, gain(count,:))
    ylabel('g(\phi)/dB')
    xlabel('\phi')
    title('Angular response versus angle');
    axis([-90 90 -350 0]);
    grid on
    hold on
end
% legend('M=8,b=12','M=8,b=4','M=32,b=12','M=32,b=4')
legend('M=8,\beta_3=0','M=8,\beta_3=-133','M=32,\beta_3=0','M=32,\beta_3=-133')

function y = tr_laprnd(angle, N)
    y = zeros(N,1);
    for iii = 1:N
        y(iii) = exp(-j*2*pi*(iii-1)*0.5*sin(angle));
    end
end