function [W_R,W_B] = MMSE_omp(W_MMSE ,A_r ,CovRx ,Mr )
%OMP�㷨���ն���С���������ջ�
W_R = [];
W_res = W_MMSE;
    for i=1:Mr
        Psi = A_r'*CovRx*W_res;      %������ؾ���
        D = diag(Psi*Psi');
        [~, k] = max(D);       %���ƥ��Ϊ��k��
        W_R = [W_R A_r(:,k)];
        W_B = (W_R'*CovRx*W_R)\W_R'*CovRx*W_MMSE;
        W_res = (W_MMSE-W_R*W_B)/norm(W_MMSE-W_R*W_B,'fro');
    end
end

