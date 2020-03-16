%__________________________________________________________________________
% violin_h.m - Simple violin plot using matlab default kernel density estimation
%
% Updates:
% v2: extended for accepting also cells of different length
% v3: 
%    - changed varargin to parameter - value list
%    - specification of x-axis vals possible now
% Last update: 09/2014
%__________________________________________________________________________
% This function creates violin plots based on kernel density estimation
% using ksdensity with default settings. Please be careful when comparing pdfs 
% estimated with different bandwidth!
%
% Differently to other boxplot functions, you may specify the x-position.
% This is particularly usefule when overlaying with other data / plots.
%__________________________________________________________________________
% 
% Please cite this function as:
% Hoffmann H, 2013: violin.m - Simple violin plot using matlab default kernel 
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
%
%__________________________________________________________________________
%
% Input:
% Y:     Data to be plotted, being either
% n x m matrix. A 'violin' is plotted for each column m, OR
% 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor=[1 0.5 0]%FaceColor: Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor='k'      %LineColor: Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha=0.5     %Alpha value (transparency)   
% mc='k'      %Color of the bars indicating the mean; set either [],'' or 'none' if the mean should not be plotted
% medc='r'    %Color of the bars indicating the median; set either [],'' or 'none' if the mean should not be plotted
% bw=[];      %Kernel bandwidth, prescribe if wanted. 
%            %If b is a single number, b will be applied to all estimates
%            %If b is an array of 1xm or mx1, b(i) will be applied to
%            column (i).
%
% Output:
% h: figure handle
% L: Legend handle
% MX: Means of groups
% MED: Medians of groups
% bw: bandwidth of kernel
%__________________________________________________________________________
%{
% Example1 (default):

disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
Y2=[rand(1000,1)+1,gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)]+10;
figure;
[h,L,MX,MED]=violin_h(Y,Y2); 
figure;
[h,L,MX,MED]=violin(Y); 
figure;
[h,L,MX,MED]=violin(Y2); 
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%Example2 (specify facecolor, edgecolor, xlabel):

disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
Y2=[rand(1000,1)+1,gamrnd(1,2,1000,1)+1,normrnd(10,2,1000,1)+1,gamrnd(10,0.1,1000,1)+1]+10;
violin_h(Y,Y2,'xlabel',{'a','b','c','d'},'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--')
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%Example3 (specify x axis location):

disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)]+1;
violin_h(Y,Y2,'x',[-1 .7 3.4 8.8],'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none',...
'bw',0.3,'mc','k','medc','r-.')
axis([-2 10 -0.5 20])
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%Example4 (Give data as cells with different n):

disp('this example uses the statistical toolbox')

Y{:,1}=rand(10,1);
Y{:,2}=rand(1000,1);
Y2{:,1}=rand(10,1)+1;
Y2{:,2}=rand(1000,1)+1;
violin_h(Y,Y2,'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none','bw',0.1,'mc','k','medc','r-.')
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%}
%__________________________________________________________________________
%__________________________________________________________________________

function[h,L,MX,MED,bw]=violin_h(Y1,Y2,varargin)

 %defaults:
 %_____________________
 xL=[];
 fc=[1 0.5 0];
 fc=[0.3 0.3 0.3];
 lc='k';
 alp=1;
 mc='k';
 medc='r';
 b=[]; %bandwidth
 b2=[]; %bandwidth
 plotlegend=0; 
 plotmean=1;
 plotmedian=0;
 plotci=1;
 side='onesided';
 cip=0.95;
 x = [];
 %_____________________

 %convert single columns to cells:
 if iscell(Y1)==0
    Y1 = num2cell(Y1,1);
 end
 if iscell(Y2)==0
    Y2 = num2cell(Y2,1);
 end
 
 %get additional parameters
 if isempty(find(strcmp(varargin,'xlabel')))==0
     xL = varargin{find(strcmp(varargin,'xlabel'))+1};
 end
 if isempty(find(strcmp(varargin,'facecolor')))==0
   fc = varargin{find(strcmp(varargin,'facecolor'))+1};
 end
 if isempty(find(strcmp(varargin,'edgecolor')))==0
     lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
 end
 if isempty(find(strcmp(varargin,'facealpha')))==0
     alp = varargin{find(strcmp(varargin,'facealpha'))+1};
 end
 if isempty(find(strcmp(varargin,'mc')))==0
   if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
     mc = varargin{find(strcmp(varargin,'mc'))+1};
     plotmean = 1;
   else
     plotmean = 0;
   end
 end
 if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
     medc = varargin{find(strcmp(varargin,'medc'))+1};
     plotmedian = 1;
    else
     plotmedian = 0;
    end
 end
 if isempty(find(strcmp(varargin,'side')))==0
     side = varargin{find(strcmp(varargin,'side'))+1};
 end
 if isempty(find(strcmp(varargin,'ci')))==0
     plotci = varargin{find(strcmp(varargin,'ci'))+1};
 end
 if isempty(find(strcmp(varargin,'cip')))==0
     cip = varargin{find(strcmp(varargin,'cip'))+1};
 end
 if isempty(find(strcmp(varargin,'bw')))==0
     b = varargin{find(strcmp(varargin,'bw'))+1}
      if length(b)==1
         disp(['same bandwidth bw = ',num2str(b),' used for all cols'])
         b=repmat(b,size(Y1,2),1);
      elseif length(b)~=size(Y1,2)
         warning('length(b)~=size(Y1,2)')
         error('please provide only one bandwidth or an array of b with same length as columns in the data set')
      end
     b2 = varargin{find(strcmp(varargin,'bw'))+1}
      if length(b2)==1
         disp(['same bandwidth bw = ',num2str(b2),' used for all cols'])
         b2=repmat(b2,size(Y2,2),1);
      elseif length(b2)~=size(Y2,2)
         warning('length(b)~=size(Y2,2)')
         error('please provide only one bandwidth or an array of b with same length as columns in the data set')
      end
 end
 if isempty(find(strcmp(varargin,'plotlegend')))==0
     plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
 end
 if isempty(find(strcmp(varargin,'x')))==0
     x = varargin{find(strcmp(varargin,'x'))+1};
 end
 %-------------------------------------------------------------------------
 if size(fc,1)==1
     fc=repmat(fc,size(Y1,2),1);
     fc2=repmat(1-fc,size(Y2,2),1);
 end
 %-------------------------------------------------------------------------
 i=1;
 for i=1:size(Y1,2)
     
     if isempty(b)==0
         [f, u, bb]=ksdensity(Y1{i},'bandwidth',b(i));
     elseif isempty(b)
         [f, u, bb]=ksdensity(Y1{i});
     end
     
     f=f/max(f)*0.3; %normalize
     F(:,i)=f;
     U(:,i)=u;
     MED(:,i)=nanmedian(Y1{i});
     MX(:,i)=nanmean(Y1{i});
     bw(:,i)=bb;
 end
 for i=1:size(Y2,2)
     
     if isempty(b2)==0
         [f2, u2, bb2]=ksdensity(Y2{i},'bandwidth',b2(i));
     elseif isempty(b2)
         [f2, u2, bb2]=ksdensity(Y2{i});
     end
     
     f2=f2/max(f2)*0.3; %normalize
     F2(:,i)=f2;
     U2(:,i)=u2;
     MED2(:,i)=nanmedian(Y2{i});
     MX2(:,i)=nanmean(Y2{i});
     bw2(:,i)=bb2;
     
 end
 %-------------------------------------------------------------------------
 %mp = get(0, 'MonitorPositions');
 %set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
 %-------------------------------------------------------------------------
 if isempty(x)
     x = zeros(size(Y1,2));
     setX = 0;
 else
  setX = 1;
  if isempty(xL)==0
   disp('_________________________________________________________________')
   warning('Function is not designed for x-axis specification with string label')
   warning('when providing x, xlabel can be set later anyway')
   error('please provide either x or xlabel. not both.')
  end
 end
 %_________________________________________________________________________
% pos=[    0.8571    1.1429]-1;
 %_________________________________________________________________________
 i=1;
 for i=i:size(Y1,2)
     if plotci==1
%          ci1=ciestim3(Y1{i},1,cip,side);
%          ci2=ciestim3(Y2{i},1,cip,side);
         if size(Y1{i}(~isnan(Y1{i})),1)==1
         ci1=ciestim3(Y1{i}(~isnan(Y1{i})),2,cip,side)';
         ci2=ciestim3(Y2{i}(~isnan(Y2{i})),2,cip,side)';
         else
         ci1=ciestim3(Y1{i}(~isnan(Y1{i})),1,cip,side);
         ci2=ciestim3(Y2{i}(~isnan(Y2{i})),1,cip,side);
         end
     end
   if isempty(lc) == 1
    if setX == 0
%      h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
     h(i)=fill([F(:,i)+i],[U(:,i)],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
     hold on
     h(i)=fill([flipud(i-F2(:,i))],[flipud(U2(:,i))],fc2(i,:),'FaceAlpha',alp,'EdgeColor','none');
    else
%      h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
     h(i)=fill([F(:,i)+x(i)],[U(:,i)],fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
     hold on
     h(i)=fill([flipud(x(i)-F2(:,i))],[flipud(U2(:,i))],fc2(i,:),'FaceAlpha',alp,'EdgeColor','none');
    end
   else
    if setX == 0
     %h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
     %h(i)=fill([F2(:,i)+i;flipud(i-F2(:,i))],[U2(:,i);flipud(U2(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
     h(i)=fill([flipud(i-F(:,i))],[flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
     hold on
     h(i)=fill([F2(:,i)+i],[U2(:,i)],fc2(i,:),'FaceAlpha',alp,'EdgeColor',lc);
    else
%      h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
     h(i)=fill([flipud(x(i)-F(:,i))],[flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
     hold on
     h(i)=fill([F2(:,i)+x(i)],[U2(:,i)],fc2(i,:),'FaceAlpha',alp,'EdgeColor',lc);
    end
     if plotci==1
     plot([i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)+ci1) ],[MX(:,i) MX(:,i)]+ci1,'r','LineWidth',1);
     plot([i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)-ci1) ],[MX(:,i) MX(:,i)]-ci1,'r','LineWidth',1);
     plot([interp1(U2(:,i),F2(:,i)+i,MX2(:,i)+ci2), i ],[MX2(:,i) MX2(:,i)]+ci2,'r','LineWidth',1);
     plot([interp1(U2(:,i),F2(:,i)+i,MX2(:,i)-ci2), i ],[MX2(:,i) MX2(:,i)]-ci2,'r','LineWidth',1);
     end
   end
    hold on
   if setX == 0
    if plotmean == 1
%      p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',1);
     p(1)=plot([i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',1);
     hold on
     p(1)=plot([interp1(U2(:,i),F2(:,i)+i,MX2(:,i)), i ],[MX2(:,i) MX2(:,i)],mc,'LineWidth',1);
    end
    if plotmedian == 1
%      p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i)) ],[MED(:,i) MED(:,i)],medc,'LineWidth',1);
     p(2)=plot([i interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i)) ],[MED(:,i) MED(:,i)],medc,'LineWidth',1);
     hold on
     p(2)=plot([interp1(U2(:,i),F2(:,i)+i,MED2(:,i)), i],[MED2(:,i) MED2(:,i)],medc,'LineWidth',1);
    end
   elseif setX == 1
    if plotmean == 1
%      p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',1);
     p(1)=plot([i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',1);
     hold on
     p(1)=plot([interp1(U2(:,i),F2(:,i)+i,MX2(:,i))+x(i)-i, i],[MX2(:,i) MX2(:,i)],mc,'LineWidth',1);
    end
    if plotmedian == 1
%      p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',1);
     p(2)=plot([i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',1);
     hold on
     p(2)=plot([interp1(U2(:,i),F2(:,i)+i,MED2(:,i))+x(i)-i, i],[MED2(:,i) MED2(:,i)],medc,'LineWidth',1);
    end
   end
 end
 %-------------------------------------------------------------------------
 if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
  
  if plotmean==1 & plotmedian==1
   L=legend([p(1) p(2)],'Mean','Median');
  elseif plotmean==0 & plotmedian==1
   L=legend([p(2)],'Median');
  elseif plotmean==1 & plotmedian==0
   L=legend([p(1)],'Mean');
  end 
  
  set(L,'box','off','FontSize',14)
 else
    L=[];
 end
 %-------------------------------------------------------------------------
 axis([0.5 size(Y1,2)+0.5, min([U(:);U2(:)]) max([U(:);U2(:)])]);
 %-------------------------------------------------------------------------
 xL2={''};
 i=1;
 for i=1:size(xL,2)
  xL2=[xL2,xL{i},{''}];
 end
 set(gca,'TickLength',[0 0],'FontSize',12)
 box on
 
 if isempty(xL)==0
     set(gca,'XtickLabel',xL2)
 end
 %-------------------------------------------------------------------------
end %of function