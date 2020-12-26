function [F_R, F_B] = omp(F_opt, A_t, Ns, Mt)
F_R = [];
F_res = F_opt;
    for i=1:Mt
        Psi = A_t'*F_res;      %计算相关矩阵
        D = diag(Psi*Psi');
        [~, k] = max(D);       %最佳匹配为第k列
        F_R = [F_R A_t(:,k)];
        F_B = (F_R'*F_R)\F_R'*F_opt;
        F_res = (F_opt-F_R*F_B)/norm(F_opt-F_R*F_B,'fro');
    end
F_B = sqrt(Ns)*F_B/norm(F_R*F_B,'fro');
end
