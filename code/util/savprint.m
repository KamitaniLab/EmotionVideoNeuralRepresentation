function verb=savprint(handle,savName,verbose)
% function savprint(handle,savName,verbose)
% 
% 
% Tomoyasu Horikawa horikawa-t@atr.jp 2010/10/26
% 


setdir(fileparts(savName));
set(handle, 'PaperType', 'A4');
set(handle,'PaperOrientation','landscape');
set(handle,'PaperUnits','normalized');
set(handle,'PaperPosition',[0 0 1 1]);
if strcmp(savName(end-2:end),'png')
    print(handle,'-dpng',savName);
end
if strcmp(savName(end-2:end),'pdf')
%     print(handle,'-dpdf',savName);
    print(handle,'-painters','-dpdf',savName);
%     print(handle,'-painters','-r300','-dpdf',savName);
end
if strcmp(savName(end-3:end),'tiff')
    print(handle,'-dtiff',savName);
end
if strcmp(savName(end-1:end),'ps')
    print(handle,'-dps',savName);
end
if strcmp(savName(end-3:end),'jpeg')||strcmp(savName(end-3:end),'jpg')
    print(handle,'-djpeg',savName);
end
verb=sprintf(savName);
if exist('verbose','var') && verbose
    fprintf(verb)
    fprintf('\n')
end
%%
% hf = gcf;
% hf.Renderer = 'painters';
% print('im.pdf','-dpdf','-bestfit')