function hh=shiftebar_h(mean,sd,nshift,varargin)
% shiftebar -- draw shifted errorbar_h
%function shiftebar(mean,sd,nshift,varargin)
%  
%  [Input]
%     -mean: mean of the data
%     -sd: std of the data
%     -nshift: number of shift
%  
%  
%  
%  [note]
% 2 [-0.14, 0.14]
% 3 [-0.14, 0, 0.14]
% 4 [-0.275, -0.095, 0.09, 0.275]
% 5 [-0.305, -0.1525, 0, 0.1525, 0.305]
% 6 [-0.33, -0.205,-0.07, 0.07, 0.205, 0.33]
%
%
%
% Created By Tomoyasu Horikawa horikawa-t@atr.jp 2010/04/19
%
%
%
% n=18;
% hh=bar(rand(n));
% xxx=zeros(n,1);
% a=get(hh,'Children');for i=1:n,b=get(a{i});xxx(i)=mean(b.XData(2:3,1));end
% xxx'

pos=cell(10,1);
pos{1}=[     1]-1;
pos{2}=[    0.8571    1.1429]-1;
pos{3}=[    0.7778    1.0000    1.2222]-1;
pos{4}=[    0.7273    0.9091    1.0909    1.2727]-1;
pos{5}=[    0.6923    0.8462    1.0000    1.1538    1.3077]-1;
pos{6}=[    0.6667    0.8000    0.9333    1.0667    1.2000    1.3333]-1;
pos{7}=[    0.6571    0.7714    0.8857    1.0000    1.1143    1.2286    1.3429]-1;
pos{8}=[    0.6500    0.7500    0.8500    0.9500    1.0500    1.1500    1.2500    1.3500]-1;
pos{9}=[    0.6444    0.7333    0.8222    0.9111    1.0000    1.0889    1.1778    1.2667    1.3556]-1;
pos{10}=[    0.6400    0.7200    0.8000    0.8800    0.9600    1.0400    1.1200    1.2000    1.2800    1.3600]-1;
pos{11}=[    0.6364    0.7091    0.7818    0.8545    0.9273    1.0000    1.0727    1.1455    1.2182    1.2909    1.3636]-1;
pos{12}=[    0.6333    0.7000    0.7667    0.8333    0.9000    0.9667    1.0333    1.1000    1.1667    1.2333    1.3000    1.3667]-1;
pos{13}=[    0.6308    0.6923    0.7538    0.8154    0.8769    0.9385    1.0000    1.0615    1.1231    1.1846    1.2462    1.3077    1.3692]-1;
pos{14}=[    0.6286    0.6857    0.7429    0.8000    0.8571    0.9143    0.9714    1.0286    1.0857    1.1429    1.2000    1.2571    1.3143    1.3714]-1;
pos{15}=[    0.6267    0.6800    0.7333    0.7867    0.8400    0.8933    0.9467    1.0000    1.0533    1.1067    1.1600    1.2133    1.2667    1.3200    1.3733]-1;
if nshift>15 
    error('Not yet implemented')
end
hh=cell(nshift,1);
for i=1:nshift
hh{i}=errorbar_h((1:size(mean,1))+pos{nshift}(i),mean(:,i),sd(:,i),varargin{:});hold on
end
%{
switch nshift
    case 1
        hh=errorbar_h((1:size(mean,1)),mean(:,1),sd(:,1),varargin{:});
    case 2
        hh{1}=errorbar_h((1:size(mean,1))-0.144,mean(:,1),sd(:,1),varargin{:});hold on
        hh{2}=errorbar_h((1:size(mean,1))+0.144,mean(:,2),sd(:,2),varargin{:});
    case 3
        hh{1}=errorbar_h((1:size(mean,1))-0.225,mean(:,1),sd(:,1),varargin{:});hold on
        hh{2}=errorbar_h((1:size(mean,1))+0,mean(:,2),sd(:,2),varargin{:});
        hh{3}=errorbar_h((1:size(mean,1))+0.225,mean(:,3),sd(:,3),varargin{:});
    case 4
        hh{1}=errorbar_h((1:size(mean,1))-0.275,mean(:,1),sd(:,1),varargin{:});hold on
        hh{2}=errorbar_h((1:size(mean,1))-0.0925,mean(:,2),sd(:,2),varargin{:}); 
        hh{3}=errorbar_h((1:size(mean,1))+0.0925,mean(:,3),sd(:,3),varargin{:}); 
        hh{4}=errorbar_h((1:size(mean,1))+0.275,mean(:,4),sd(:,4),varargin{:});
    case 5
        hh{1}=errorbar_h((1:size(mean,1))-0.305,mean(:,1),sd(:,1),varargin{:});hold on
        hh{2}=errorbar_h((1:size(mean,1))-0.1525,mean(:,2),sd(:,2),varargin{:});
        hh{3}=errorbar_h((1:size(mean,1))+0,mean(:,3),sd(:,3),varargin{:});
        hh{4}=errorbar_h((1:size(mean,1))+0.1525,mean(:,4),sd(:,4),varargin{:});
        hh{5}=errorbar_h((1:size(mean,1))+0.305,mean(:,5),sd(:,5),varargin{:});
    case 6
        hh{1}=errorbar_h((1:size(mean,1))-0.33,mean(:,1),sd(:,1),varargin{:});hold on
        hh{2}=errorbar_h((1:size(mean,1))-0.205,mean(:,2),sd(:,2),varargin{:});
        hh{3}=errorbar_h((1:size(mean,1))-0.07,mean(:,3),sd(:,3),varargin{:});
        hh{4}=errorbar_h((1:size(mean,1))+0.07,mean(:,4),sd(:,4),varargin{:});
        hh{5}=errorbar_h((1:size(mean,1))+0.205,mean(:,5),sd(:,5),varargin{:});
        hh{6}=errorbar_h((1:size(mean,1))+0.33,mean(:,6),sd(:,6),varargin{:});
    case 7
        nshift=7;
        dif=0.114;
        len=size(mean,1);
        for i=1:nshift
            hh{i}=errorbar_h((1:len)'-(dif*nshift/2-i*dif)-0.055,mean(:,i),sd(:,i),varargin{:});hold on
            set(hh{i},'Color',[0,0,0]);
        end
    case 8
        nshift=8;
        dif=0.1;
        len=size(mean,1);
        for i=1:nshift
            hh{i}=errorbar_h((1:len)'-(dif*nshift/2-i*dif)-0.049,mean(:,i),sd(:,i),varargin{:});hold on
            set(hh{i},'Color',[0,0,0]);
        end
    case 9
        nshift=9;
        dif=0.089;
        len=size(mean,1);
        for i=1:nshift
            hh{i}=errorbar_h((1:len)'-(dif*nshift/2-i*dif)-0.042,mean(:,i),sd(:,i),varargin{:});hold on
            set(hh{i},'Color',[0,0,0]);
        end
    case 10
        nshift=10;
        dif=0.08;
        len=size(mean,1);
        for i=1:nshift
            hh{i}=errorbar_h((1:len)'-(dif*nshift/2-i*dif)-0.04,mean(:,i),sd(:,i),varargin{:});hold on
            set(hh{i},'Color',[0,0,0]);
        end
    case 11
        nshift=11;
        dif=0.073;
        len=size(mean,1);
        for i=1:nshift
            hh{i}=errorbar_h((1:len)'-(dif*nshift/2-i*dif)-0.037,mean(:,i),sd(:,i),varargin{:});hold on
            set(hh{i},'Color',[0,0,0]);
        end
    case 13
        nshift=13;
        dif=0.062;
        len=size(mean,1);
        for i=1:nshift
            hh{i}=errorbar_h((1:len)'-(dif*nshift/2-i*dif)-0.03,mean(:,i),sd(:,i),varargin{:});hold on
            set(hh{i},'Color',[0,0,0]);
        end
    otherwise
        hh=errorbar_h(mean);
end
%}
if nshift >1
hh=merge(hh);
end

%%
