function [Temp_ber] = ber_4QAM(S, Ns, data1)
%接收端误码率计算
S_re=real(S);
S_im=imag(S);
Temp_ber_begin=0;
for k=1:Ns
    if S_re(k,1)>0
       S_re(k,1)=1;
    else
       S_re(k,1)=-1;
    end
    if S_im(k,1)>0
       S_im(k,1)=1;
    else
       S_im(k,1)=-1;
    end
    if S_re(k,1)~=data1(k,1)
        Temp_ber_begin=Temp_ber_begin+1;
    end
    if S_im(k,1)~=data1(k,2)
        Temp_ber_begin=Temp_ber_begin+1;
    end
end
Temp_ber = Temp_ber_begin;
end

