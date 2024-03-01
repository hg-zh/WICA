function [out]= PUT_F(x,j,k)
    if (x<=1&x>=0)
        PUT_q = evalin('base', 'PUT_mean_q');
        out=interp1([0:0.01:1],PUT_q(j,:),x)-(k-1)/100;
    elseif (x>1)
        out=1;
    else
        out=-1;
    end
end