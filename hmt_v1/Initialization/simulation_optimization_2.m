
simnames={'well_2_2'};
thresh=.12; %increase if nucleus not included
dil=2; %increase if nucleus or disconnected

meth = csvread([simnames{1},'_me3.csv'],1,0);
acet = csvread([simnames{1},'_ac.csv'],1,0);
px=10;
xm=meth(:,3); ym=meth(:,4); zm=meth(:,5);
xa=acet(:,3); ya=acet(:,4); za=acet(:,5);
minx=min([xm;xa]); miny=min([ym;ya]); minz=min([zm;za]);
xm2=(xm+.01-minx)/px; ym2=(ym+.01-miny)/px; zm2=(zm+.01-(minz))/px;%(pixel size will be 8nm; you can choose other number)
xa2=(xa+.01-minx)/px; ya2=(ya+.01-miny)/px; za2=(za+.01-(minz))/px;%(pixel size will be 8nm; you can choose other number)

xcomb=[xm2;xa2]; ycomb=[ym2;ya2]; zcomb=[zm2;za2];
a1=length(xm2);
imeth=zeros(floor(max(xcomb)+5),floor(max(ycomb)+5),floor(max(zcomb)+5));

for ind=1:a1
   imeth(floor(xm2(ind))+1,floor(ym2(ind))+1,floor(zm2(ind))+1)=imeth(floor(xm2(ind))+1,floor(ym2(ind))+1,floor(zm2(ind))+1)+1;
end

%Generate the image of the second split
imac=zeros(floor(max(xcomb)+5),floor(max(ycomb)+5),floor(max(zcomb)+5));
ac_hmap=imac;
imac_thresh=repmat((imac(:,:,1)),[1,1,3]);
a2=length(xa2);
for ind=1:a2
   imac(floor(xa2(ind))+1,floor(ya2(ind))+1,floor(za2(ind))+1)=imac(floor(xa2(ind))+1,floor(ya2(ind))+1,floor(za2(ind))+1)+1;

end
%load('gt_values')

cell_img2bin = imbinarize(imgaussfilt(sum(imeth+imac,3),1), thresh);
cell_img2bin = imgaussfilt(double(cell_img2bin),dil);
cell_img2bin=imfill(cell_img2bin,'holes');
cell_img2bin2 = bwareaopen(logical(cell_img2bin),(1000000));
cell_img2bin=cell_img2bin2.*cell_img2bin;
figure; imagesc(logical(cell_img2bin))
sa=sum(sum(logical(cell_img2bin)))*(.01^2);
%%
ztot=[zm;za];h=histogram(ztot,'BinWidth',10);
mid=find(h.Values==max(h.Values));
lrng=h.BinEdges(mid)-500;
hrng=h.BinEdges(mid)+500;
ztot(ztot<lrng)=[]; ztot(ztot>hrng)=[]; h=histogram(ztot,'BinWidth',10);
N=h.Values/sum(h.Values);
for n=1:length(N)
frame(n)=round((h.BinEdges(n)+h.BinEdges(n+1))/2);
end
px=10; %pixel size
nmeth=length(xm);
nac=length(xa);
load('initial_variables')
%%
[cell_img2bin,sa]=find_area(px,simnames,thresh,dil,cell_img2bin);

%%
[meth_withinc,mdistz,ac_within_ac_c,adistz]=find_within2(xm,ym,zm,xa,ya,za,px,cell_img2bin,xm2,ym2,xa2,ya2);
nseedm=length(meth_withinc);
nseeda=length(ac_within_ac_c);
%%
sim_histones_3d('Simulation',nmeth,nac,ratiom_i,ratioa_i,ratiomz_i,ratioaz_i,nseedm,nseeda,sa,N,frame,meth_betweeni,ac_between_aci)
split_me_ac
meth = csvread('Simulation_me3.csv',1,0);
acet = csvread('Simulation_ac.csv',1,0);
xm_sim=meth(:,3); ym_sim=meth(:,4); zm_sim=meth(:,5);
xa_sim=acet(:,3); ya_sim=acet(:,4); za_sim=acet(:,5);
[sim_meth_withinc,sim_mdistz,sim_ac_within_ac_c,sim_adistz]=find_within(xm_sim,ym_sim,zm_sim,xa_sim,ya_sim,za_sim,px);
[ratiom,ratioa,ratiomz,ratioaz,meth_between,ac_between_ac,dif]=create_variables(meth_withinc,mdistz,ac_within_ac_c,adistz,sim_meth_withinc,sim_mdistz,sim_ac_within_ac_c,sim_adistz,ratiom_i,ratioa_i,ratiomz_i,ratioaz_i);
iterationN=0;
%%
while (iterationN<10)&&(dif>0.3)
iterationN=iterationN+1;
sim_histones_3d('Simulation',nmeth,nac,ratiom,ratioa,ratiomz,ratioaz,nseedm,nseeda,sa,N,frame,meth_between,ac_between_ac)
split_me_ac
meth = csvread('Simulation_me3.csv',1,0);
acet = csvread('Simulation_ac.csv',1,0);
xm_sim=meth(:,3); ym_sim=meth(:,4); zm_sim=meth(:,5);
xa_sim=acet(:,3); ya_sim=acet(:,4); za_sim=acet(:,5);
[sim_meth_withinc,sim_mdistz,sim_ac_within_ac_c,sim_adistz]=find_within(xm_sim,ym_sim,zm_sim,xa_sim,ya_sim,za_sim,px);
[ratiom,ratioa,ratiomz,ratioaz,meth_between,ac_between_ac,dif]=create_variables(meth_withinc,mdistz,ac_within_ac_c,adistz,sim_meth_withinc,sim_mdistz,sim_ac_within_ac_c,sim_adistz,ratiom,ratioa,ratiomz,ratioaz);

end
%%
%%
 [ac_width,m_width,dens,densa]=find_gt_width(meth_between,nmeth,ac_between_ac,nac);
 save('sim_param','nmeth','nac','ratiom','ratioa','ratiomz','ratioaz','nseedm','nseeda','sa','N','frame','meth_between','ac_between_ac')
 save('gt_values','ac_width','m_width','dens','densa','thresh','dil')
%%