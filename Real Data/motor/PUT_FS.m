function [out]= PUT_FS(x,j,k)
    global PUT_mean_q
    if (x<=1&x>=0)
        out=interp1([0:0.01:1],PUT_mean_q(j,:),x)-(k-1)/100;
    elseif (x>1)
        out=1;
    else
        out=-1;
    end
end