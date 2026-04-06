function [xm,ym,zm]=exclude_bg(xm,ym,zm,cell_img2bin,minx,miny,px)
xm2=(xm+.01-minx)/px; ym2=(ym+.01-miny)/px;%(pixel size will be 8nm; you can choose other number)
for n1=1:length(xm)
    if cell_img2bin(floor(xm2(n1))+1,floor(ym2(n1))+1)==0
        xm(n1)=nan;
        ym(n1)=nan;
        zm(n1)=nan;
    end
end
xm(isnan(xm))=[];
ym(isnan(ym))=[];
zm(isnan(zm))=[];
