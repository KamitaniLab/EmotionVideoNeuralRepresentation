function [Tte,W,Va]=flinreg_l2(lambda,X1X1,Xtr_1,Ttr,Xte)
% linreg_ml -- perform linear regression by maximum likelihood estimate
% [Tte,W,Va]=linreg_ml(Xtr,Ttr,Xte)
%
% [Input]
%   -Xtr: training data [N x M]
%   -Ttr: target vector [N x 1]
%   -Xte: test data  [N2 x M]
%
%
% [Output]
%   -Tte: prediction from the Xte  [N2 x 1]
%   -W: weight  [M x 1]
%   -Va: variance [1 x 1]
%
%
% [Example]
%  a=3;
%  intercept=20;
%  b=100;
%  Xtr=(1:100)';
%  Ttr=[Xtr,ones(size(Xtr,1),1)]*[a,intercept]'+randn(100,1)*sqrt(b);
%  subplot(1,3,1)
%  plot(Xtr,Ttr,'b');hold on
%  Xte=(101:200)';
%  [Tte,W,Va]=linreg_ml(Xtr,Ttr,Xte)
%  [w]=regress(Ttr,addbias(Xtr))
%  [wnn]=lsqnonneg(addbias(Xtr),Ttr)
%  plot(Xtr,[Xtr,ones(size(Xtr,1),1)]*W,':r');
%  plot(Xte,Tte,'g');hold off
%  subplot(1,3,2)
%  plot(Xtr,[Xtr,ones(size(Xtr,1),1)]*W,'-r');hold on
%  plot(Xtr,Ttr,'ok');hold off
%  subplot(1,3,3)
%  plot(Xtr,[Xtr,ones(size(Xtr,1),1)]*wnn,'-r');hold on
%  plot(Xtr,Ttr,'ok');hold off
% 
%
% [Related function]
%  
%
% Created  By: Tomoyasu Horikawa horikawa-t@atr.jp 2009/10/11
% Modified by: Tomoyasu Horikawa horikawa-t@atr.jp 2011/7/15 
%  add rcond lines
%
%

% calcurate the weights
X1X1 = X1X1 + eye(size(Xtr_1,2))*lambda;
nanIdx = all(isnan(X1X1));
W = zeros(size(X1X1,2),size(Ttr,2));

if rcond(X1X1(~nanIdx,~nanIdx))>1e-8
    W(~nanIdx,:)=X1X1(~nanIdx,~nanIdx)\(Xtr_1(:,~nanIdx)'*Ttr);%inv(XX)*Xtr_1'*Ttr;
else
    W(~nanIdx,:)=pinv(X1X1(~nanIdx,~nanIdx))*Xtr_1(:,~nanIdx)'*Ttr;
end

% predict from test data
Xte1 = [Xte,ones(size(Xte,1),1)];
Tte = Xte1(:,~nanIdx)*W(~nanIdx,:);

% calcurate the variance of the prediction
if nargout==3
Va=((Ttr-Xtr_1*W)'*(Ttr-Xtr_1*W))/size(Xtr_1,1);
end




