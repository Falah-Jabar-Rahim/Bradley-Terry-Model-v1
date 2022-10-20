function [Q_log, Lb, Ub]= BT_model(P)
% Bradley-Terry-Model
% (c) Falah Jabar (falah.jabar@lx.it.pt)

%% input:
% P : wining frequency matrix
%% Output 
% Q: BT scores (or quality scores)
% Lb: lower bound of confidence interval
% Ub: upper bound of confidence interval
 
Stimu=length(P(:,1)); %% number of conditions or stimuli
%% Cost function
f = @(x)pwcp2quality(x,P,Stimu);
lb =zeros (1,Stimu);%% lowerbound
ub = ones (1,Stimu); %% upperbound
A =  [];
b = [];
Aeq = [];
beq = [];
v0=zeros (1,Stimu);%% intial values
v0(:)=1/Stimu;%% intial values
nonlcon = @unitdisk;
options = optimset( 'Display','iter','UseParallel',1,'PlotFcn','optimplotfval');
options.Algorithm = 'interior-point';
[Q,fval,exitflag,output,lambda,grad,hessian] = fmincon(f,v0,A,b,Aeq,beq,lb,ub,nonlcon,options);

%% computes numerical Hessian
[CM ,N]= HHessian(P,Q);
Sigm=diag(CM);
Ub=zeros (1,Stimu);%% intial values
Lb=zeros (1,Stimu);%% intial values

for i=1:Stimu
   Lb(i)=log(Q(i))-(1.96*sqrt(Sigm(i)/N)/Q(i));
   Ub(i)=log(Q(i))+(1.96*sqrt(Sigm(i)/N)/Q(i));
 
end
Q_log=log(Q);
%% row2col
Q_log=Q_log'; Lb=Lb'; Ub=Ub'; Q=Q';

function [c,ceq] = unitdisk(x)
c = [];
ceq = sum(x) - 1;
end

end
