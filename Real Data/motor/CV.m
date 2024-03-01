function [out]=CV(epsilon)
    Train = evalin('base', 'Train');
    Test_CAU = evalin('base', 'Test_CAU');
    Test_PUT = evalin('base', 'Test_PUT');
    S=0;
        for fold=1:length(Train)
            C=inv(Train{fold}.Ax+epsilon*eye(50))*Train{fold}.Axy*inv(Train{fold}.Ay+epsilon*eye(50))*Train{fold}.Ayx;
            [V,D]=eig(C);
            [~,ind]=max(real(diag(D)));
            R=real(V(:,ind));
            G=inv(Train{fold}.Ay+epsilon*eye(50))*Train{fold}.Ayx*R/norm(inv(Train{fold}.Ay+epsilon*eye(50))*Train{fold}.Ayx*R);
            F_hat=zeros(284,101);G_hat=zeros(284,101);
            for k=1:50
                F_hat=F_hat+R(k)*Train{fold}.CAU.basis{k};
                G_hat=G_hat+G(k)*Train{fold}.PUT.basis{k};
            end
            A_x=zeros(length(Test_CAU{fold}),1);A_y=zeros(length(Test_PUT{fold}),1);
            for i=1:length(Test_CAU{fold})
                tmp=0;
                for j=1:284
                    tmp=tmp+sum(Test_CAU{fold}{i}(j,:).*F_hat(j,:).*[Train{fold}.CAU.pdf(j,:),0] );
                end
                A_x(i)=tmp*0.01/284; tmp=0;
                for j=1:284
                    tmp=tmp+sum(Test_PUT{fold}{i}(j,:).*G_hat(j,:).*[Train{fold}.PUT.pdf(j,:),0] );
                end
                A_y(i)=tmp*0.01/284;
            end
            cor=corrcoef(A_x,A_y);
            S=S+cor(1,2)^2;
        end
    out=-S;
end