function draw_figure4F(p,encRes)
%
% This code is for drawing figure 4F
%
%
%% settings
compPairs = {'category','dimension'};

nSbj = length(p.sbjID);

%% draw figures 4F
% summarize idnetification accuracy
iden_cr_mu = cell(length(compPairs),nSbj);
iden_cr_ci = cell(length(compPairs),nSbj);
for sbjitr = 1:nSbj
    for scoritr = 1:length(compPairs)
        fset1 = compPairs{scoritr};
        acc1 = encRes.(fset1).iden_acc{sbjitr};
        
        [ci,mu] = ciestim3(acc1,1,0.99,'onesided');
        iden_cr_mu{scoritr,sbjitr} = mu;
        iden_cr_ci{scoritr,sbjitr} = ci;
        
    end
end

% display settings
r = 5;
c = 5;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[2,1]);
fsize = 8;
msize = 2;%1.5;

close all
h = ffigure;
cmap3('ibg4');
cnt = 0;


cnt = cnt + 1;
subplottight(r,c,o(cnt),0.1);

% summary
mu = zeros(nSbj,length(compPairs));
ci = zeros(nSbj,length(compPairs));
for ix = 1:length(compPairs)
    for sbjitr = 1:nSbj
        scoreidx = ismember(compPairs,compPairs{ix});
        mu(sbjitr,ix) = iden_cr_mu{scoreidx,sbjitr}*100;
        ci(sbjitr,ix) = iden_cr_ci{scoreidx,sbjitr}*100;
    end
end

hold on
bar(mean(mu,1)','edgecolor','none')
errorbar_h(mean(mu,1)',ciestim3(mu,1,0.99,'onesided')','.k');
dif = 0.2;
xloc = repmat((1:length(compPairs))+dif,size(mu,1),1);
xloc = xloc+repmat(randn(size(xloc,1),1)/50,1,2);
plot(xloc',mu','-o','Color',[1,1,1]*0.4,'MarkerEdgeColor',[1,1,1]*0.2,'MarkerFaceColor',[1,1,1]*0.4,'MarkerSize',msize);
set(gca,'FontSize',fsize)
ylabel(sprintf('Identification accuracy (%%)'))
axname(compPairs,1)
ylim([40,100])
xlim([0,6])
xticklabel_rotate([],45)
hline(50,'-k')
hh = hline(40:10:100,'-k');
set(hh,'Color',[1,1,1]*0.8)
axis square
ffine(h)

suptitle(sprintf('Mean identification accuracy via encoding analyses'));
fname = [p.figdir,'figure4F.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));




%%
close all
