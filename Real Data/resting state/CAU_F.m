function [out]= CAU_F(x,j,k)
    if (x<=1&x>=0)
        CAU_q = evalin('base', 'CAU_mean_q');
        out=interp1([0:0.01:1],CAU_q(j,:),x)-(k-1)/100;
    elseif (x>1)
        out=1;
    else
        out=-1;
    end
end