function AnnotationOff(hh,idx)
for i=1:length(idx)
    x=get(hh(idx(i)),'Annotation');
    x.LegendInformation.IconDisplayStyle='off';
end
