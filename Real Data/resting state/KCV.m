function [out]=KCV(M)
        Train = evalin('base', 'Train');
        Test_CAU = evalin('base', 'Test_CAU');
        Test_PUT = evalin('base', 'Test_PUT');
        S=zeros(M,1);
        for fold=1:length(Train)
            for m=1:M
                Ax=Train{fold}.Ax(1:m,1:m);Ay=Train{fold}.Ay(1:m,1:m);
                Axy=Train{fold}.Axy(1:m,1:m);Ayx=Train{fold}.Ayx(1:m,1:m);
                C=inv(Ax)*Axy*inv(Ay)*Ayx;
                [V,D]=eig(C);
                [~,ind]=max(real(diag(D)));
                R=real(V(:,ind));
                G=inv(Ay)*Ayx*R/norm(inv(Ay)*Ayx*R);
                F_hat=zeros(284,101);G_hat=zeros(284,101);
                for k=1:m
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
                S(m)=S(m)+cor(1,2)^2;
            end
        end
        [~,out]=max(S);
end