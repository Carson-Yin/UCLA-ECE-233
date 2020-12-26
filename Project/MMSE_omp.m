function [W_R,W_B] = MMSE_omp(W_MMSE ,A_r ,CovRx ,Mr )
%OMP算法接收端最小均方误差接收机
W_R = [];
W_res = W_MMSE;
    for i=1:Mr
        Psi = A_r'*CovRx*W_res;      %计算相关矩阵
        D = diag(Psi*Psi');
        [~, k] = max(D);       %最佳匹配为第k列
        W_R = [W_R A_r(:,k)];
        W_B = (W_R'*CovRx*W_R)\W_R'*CovRx*W_MMSE;
        W_res = (W_MMSE-W_R*W_B)/norm(W_MMSE-W_R*W_B,'fro');
    end
end

