simnames={'simdens1surround2','simdens2surround2','simdens3surround2','simdens4surround2','simdens5surround2','simdens6surround2','simdens7surround2','simdens8surround2','simdens9surround2','simdens10surround2'};
n=1;
meth = csvread([simnames{n},'_me3.csv'],1,0);
acet = csvread([simnames{n},'_ac.csv'],1,0);
px=10;
xm=meth(:,3); ym=meth(:,4); zm=meth(:,5);
xa=acet(:,3); ya=acet(:,4); za=acet(:,5);
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
load('gt_values')

cell_img2bin = imbinarize(imgaussfilt(sum(imeth+imac,3),1), .17);
cell_img2bin = imgaussfilt(double(cell_img2bin),2);
cell_img2bin=imfill(cell_img2bin,'holes');
cell_img2bin2 = bwareaopen(logical(cell_img2bin),(100000));
cell_img2bin=cell_img2bin2.*cell_img2bin;
figure; imagesc(logical(cell_img2bin))
sa=sum(sum(logical(cell_img2bin)))*(.01^2);
%%
for n1=1:length(xm)
    if cell_img2bin(floor(xm(n1))+1,floor(ym(n1))+1)==0
        xm(n1)=nan;
        ym(n1)=nan;
    end
end
for n2=1:length(xa)
    if cell_img2bin(floor(xa(n2))+1,floor(ya(n2))+1)==0
        xa(n2)=nan;
        ya(n2)=nan;
    end
end
xm(isnan(xm))=[];
ym(isnan(ym))=[];
xa(isnan(xa))=[];
ya(isnan(ya))=[];
%%
epsilon=[];
densitym=[];
gt=m_width;
factor=32;
for n=[1,2,5,9,10]
meth = csvread([simnames{n},'_me3.csv'],1,0);
acet = csvread([simnames{n},'_ac.csv'],1,0);
px=10;
xm=meth(:,3); ym=meth(:,4); zm=meth(:,5);
%[xm,ym,zm]=exclude_bg(xm,ym,zm,cell_img2bin,minx,miny,px);
factor0=30;
iterationN=1;
factor=factor-2;
d_change=0;
dir=1;
cost=100;
while (abs(cost)>1.3)&&(iterationN<20)
    
    cost = epsilon_cost(factor,xm,ym,zm,gt);
    if (cost<0)&&(abs(cost)>1.3)
        if (dir==-1)&&(iterationN>1)
           d_change=d_change+1; 
        end
        factor=factor+5/(2^d_change);
        dir=1;
    elseif (cost>0)&&(abs(cost)>1.3)
        if (dir==1)&&(iterationN>1)
           d_change=d_change+1; 
        end
        factor=factor-5/(2^d_change);
        dir=-1;
    end
    iterationN=iterationN+1;
end
epsilon=[epsilon;factor]; 
densitym=[densitym;length(xm)/sa];
end
par=polyfit(log10(densitym),log10(epsilon),1);
figure; plot(densitym,epsilon,'.');
exp=par(1);
coef=10^(par(2));
save('methyl_epsilon_surround2','coef','exp');
%%
epsilona=[];
densitya=[];
gta=ac_width;
factor=32;
for n=[1,2,5,9,10]
acet = csvread([simnames{n},'_ac.csv'],1,0);
px=10;
xa=acet(:,3); ya=acet(:,4); za=acet(:,5);
%[xa,ya,za]=exclude_bg(xa,ya,za,cell_img2bin,minx,miny,px);
factor0=30;
iterationN=1;
factor=factor-2;
d_change=0;
dir=1;
cost=100;
while (abs(cost)>1.3)&&(iterationN<20)
    
    cost = epsilon_cost(factor,xa,ya,za,gta);
    if (cost<0)&&(abs(cost)>1.3)
        if (dir==-1)&&(iterationN>1)
           d_change=d_change+1; 
        end
        factor=factor+5/(2^d_change);
        dir=1;
    elseif (cost>0)&&(abs(cost)>1.3)
        if (dir==1)&&(iterationN>1)
           d_change=d_change+1; 
        end
        factor=factor-5/(2^d_change);
        dir=-1;
    end
    iterationN=iterationN+1;
end
epsilona=[epsilona;factor]; 
densitya=[densitya;length(xa)/sa];
end
%%
for n=1:10
meth = csvread([simnames{n},'_me3.csv'],1,0);
acet = csvread([simnames{n},'_ac.csv'],1,0);
xa=acet(:,3); ya=acet(:,4); za=acet(:,5);
xm=meth(:,3); ym=meth(:,4); zm=meth(:,5);

[decaym(n),decaya(n)]=calc_decay(xa,ya,za,xm,ym,zm);
densitya2(n)=length(xa)/sa;
densitym2(n)=length(xm)/sa;
end
%%
par=polyfit(log10(densitya),log10(epsilona),1);
figure; plot(densitya,epsilona,'.');
expa=par(1);
coefa=10^(par(2));
save('acetyl_epsilon_surround2','coefa','expa');
%%
