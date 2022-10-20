
function [C_mod_inv, V]= HHessian(P,x)  % computes numerical Hessian
% based on a solution posted on Matlab Central by Paul L. Fackler

i_row=length(P(:,1));
j_col=length(P(1,:));

N=0;
for i=1:i_row
    for j=1:j_col
        if i<j
            N= N+(P(i,j)+P(j,i));
        end
    end
end
Mu=zeros(i_row,j_col);
for i=1:i_row   
    for j=1:j_col            
    Mu(i,j)=(P(i,j)+P(j,i))/N;            
    end       
end

lam=zeros(i_row,j_col);
for i=1:i_row
    for j=1:j_col
        if i~=j
            lam(i,j)= -Mu(i,j)/((x(i)+x(j))^2);
        else
            lam(i,j)=0;
            for k=1:j_col
                if i~=k                                      
                    lam(i,j) = lam(i,j) + (Mu(i,k)*x(k))/((x(i)+x(k))^2);
                end
            end
            lam(i,j) = lam(i,j)/x(i);
        end
    end
end

C_mod = ([lam ,ones(i_row,1); ones(1,i_row) ,0]);
C_mod_inv=inv(C_mod);
V=N;
