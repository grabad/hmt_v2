function [cell_img2bin,sa]=find_area(px,title,thresh,dil,cell_img2bin)
n=10;
meth = csvread([title{1},'_me3.csv'],1,0);
acet = csvread([title{1},'_ac.csv'],1,0);
xm=meth(:,3); ym=meth(:,4); zm=meth(:,12);
xa=acet(:,3); ya=acet(:,4); za=acet(:,12);
minx=min([xm;xa]); miny=min([ym;ya]); minz=min([zm;za]);
xm=(xm+.01-minx)/px; ym=(ym+.01-miny)/px; zm=(zm+.01-(minz))/px;%(pixel size will be 8nm; you can choose other number)
xa=(xa+.01-minx)/px; ya=(ya+.01-miny)/px; za=(za+.01-(minz))/px;%(pixel size will be 8nm; you can choose other number)
xcomb=[xm;xa]; ycomb=[ym;ya]; zcomb=[zm;za];
a1=length(xm);
imeth=zeros(floor(max(xcomb)+5),floor(max(ycomb)+5),floor(max(zcomb)+5));

for ind=1:a1
   imeth(floor(xm(ind))+1,floor(ym(ind))+1,floor(zm(ind))+1)=imeth(floor(xm(ind))+1,floor(ym(ind))+1,floor(zm(ind))+1)+1;
end

%Generate the image of the second split
imac=zeros(floor(max(xcomb)+5),floor(max(ycomb)+5),floor(max(zcomb)+5));
ac_hmap=imac;
imac_thresh=repmat((imac(:,:,1)),[1,1,3]);
a2=length(xa);
for ind=1:a2
   imac(floor(xa(ind))+1,floor(ya(ind))+1,floor(za(ind))+1)=imac(floor(xa(ind))+1,floor(ya(ind))+1,floor(za(ind))+1)+1;

end

cell_img2bin = (cell_img2bin).*imbinarize(imgaussfilt(sum(imeth+imac,3),1), thresh);
cell_img2bin = imgaussfilt(double(cell_img2bin),dil);
cell_img2bin=imfill(cell_img2bin,'holes');
cell_img2bin2 = bwareaopen(logical(cell_img2bin),(100000));
cell_img2bin=cell_img2bin2.*cell_img2bin;
figure; imagesc(logical(cell_img2bin))
sa=sum(sum(logical(cell_img2bin)))*(.01^2);
end