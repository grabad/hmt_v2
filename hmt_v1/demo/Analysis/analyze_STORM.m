clear;

simnames={'k27_k27_thaw011'};

for k=1:length(simnames)
simnames={'k27_k27_thaw011'};
meth = csvread([simnames{k},'_C1.csv'],1,0);
acet = csvread([simnames{k},'_C2.csv'],1,0);

xm=meth(:,3); ym=meth(:,4); zm=meth(:,12);
[xm,order]=sort(xm); ym=ym(order); zm=zm(order);
xa=acet(:,3); ya=acet(:,4); za=acet(:,12);
[xa,order]=sort(xa); ya=ya(order); za=za(order);


meth_within=zeros(length(xm),25);
ac_within=zeros(length(xm),25);
sdis=1000;
for ind2=1:length(xm)
       mask=(abs(xm-xm(ind2))<=sdis).*(abs(ym-ym(ind2))<=sdis).*(abs(zm-zm(ind2))<=sdis);
       mask=logical(mask); 
           
       [ind,dist]=rangesearch([xm(mask),ym(mask),zm(mask)],[xm(ind2),ym(ind2),zm(ind2)],sdis);
        for nnn=10:10:sdis
            meth_within(ind2,nnn/10)=sum((dist{1}<=nnn))-1;
        end
        
        mask=(abs(xa-xm(ind2))<=sdis).*(abs(ya-ym(ind2))<=sdis).*(abs(za-zm(ind2))<=sdis);
        mask=logical(mask);    
       [ind,dist]=rangesearch([xa(mask),ya(mask),za(mask)],[xm(ind2),ym(ind2),zm(ind2)],sdis);
        for nnn=10:10:sdis
            ac_within(ind2,nnn/10)=sum((dist{1}<=nnn));
        end
end

meth_within_ac=zeros(length(xa),25);
ac_within_ac=zeros(length(xa),25);
for ind2=1:length(xa)
       mask=(abs(xm-xa(ind2))<=sdis).*(abs(ym-ya(ind2))<=sdis).*(abs(zm-za(ind2))<=sdis);
       mask=logical(mask);
  
       [ind,dist]=rangesearch([xm(mask),ym(mask),zm(mask)],[xa(ind2),ya(ind2),za(ind2)],sdis);
        for nnn=10:10:sdis
            meth_within_ac(ind2,nnn/10)=sum((dist{1}<=nnn));
        end
        
        mask=(abs(xa-xa(ind2))<=sdis).*(abs(ya-ya(ind2))<=sdis).*(abs(za-za(ind2))<=sdis);
        mask=logical(mask);
   
       [ind,dist]=rangesearch([xa(mask),ya(mask),za(mask)],[xa(ind2),ya(ind2),za(ind2)],sdis);
        for nnn=10:10:sdis
            ac_within_ac(ind2,nnn/10)=sum((dist{1}<=nnn))-1;
        end
end


sdis=1000;
pix_meth_in_meth=zeros(1,length(xm));
pix_meth_in_ac=zeros(1,length(xm));

for ind2=1:length(xm)
       mask=(abs(xm-xm(ind2))<=sdis).*(abs(ym-ym(ind2))<=sdis).*(abs(zm-zm(ind2))<=sdis);
       mask=logical(mask); 
           
       [ind,dist]=knnsearch([xm(mask),ym(mask),zm(mask)],[xm(ind2),ym(ind2),zm(ind2)],'K',10);
       pix_meth_in_meth(ind2)=mean(dist);
        
        mask=(abs(xa-xm(ind2))<=sdis).*(abs(ya-ym(ind2))<=sdis).*(abs(za-zm(ind2))<=sdis);
        mask=logical(mask);    
       [ind,dist]=knnsearch([xa(mask),ya(mask),za(mask)],[xm(ind2),ym(ind2),zm(ind2)],'K',10);
       pix_meth_in_ac(ind2)=mean(dist);
end

pix_ac_in_meth=zeros(1,length(xa));
pix_ac_in_ac=zeros(1,length(xa));
for ind2=1:length(xa)
       mask=(abs(xm-xa(ind2))<=sdis).*(abs(ym-ya(ind2))<=sdis).*(abs(zm-za(ind2))<=sdis);
       mask=logical(mask);
  
       [ind,dist]=knnsearch([xm(mask),ym(mask),zm(mask)],[xa(ind2),ya(ind2),za(ind2)],'K',10);
       pix_ac_in_meth(ind2)=mean(dist);
        
        mask=(abs(xa-xa(ind2))<=sdis).*(abs(ya-ya(ind2))<=sdis).*(abs(za-za(ind2))<=sdis);
        mask=logical(mask);
   
       [ind,dist]=knnsearch([xa(mask),ya(mask),za(mask)],[xa(ind2),ya(ind2),za(ind2)],'K',10);
        pix_ac_in_ac(ind2)=mean(dist);
end
%%
% figure; plot((10:10:sdis),mean(meth_within(:,:),1)./mean(ac_within(:,:),1),'g')
% hold on; plot((10:10:sdis),mean(ac_within_ac(:,:),1)./mean(meth_within_ac(:,:),1),'r')
% legend('methylation:acetylation contact ratio', 'acetylation:methylation contact ratio')
% xlabel('radius (nm)')
% ylabel('contact ratio')


% figure; plot((10:10:sdis),(mean(meth_within(:,:),1)./mean(ac_within(:,:),1))./(max((mean(meth_within(:,:),1)./mean(ac_within(:,:),1)))),'g')
% hold on; plot((10:10:sdis),(mean(ac_within_ac(:,:),1)./mean(meth_within_ac(:,:),1))./(max((mean(ac_within_ac(:,:),1)./mean(meth_within_ac(:,:),1)))),'r')
meth_ac_con=mean(meth_within(:,:),1)./mean(ac_within(:,:),1);
ac_meth_con=mean(ac_within_ac(:,:),1)./mean(meth_within_ac(:,:),1);
% figure; plot((10:10:sdis),(mean(meth_within(:,:),1)./mean(ac_within(:,:),1))./meth_ac_con(end),'g')
% hold on; plot((10:10:sdis),(mean(ac_within_ac(:,:),1)./mean(meth_within_ac(:,:),1))./ac_meth_con(end),'r')

% legend('Normalized methylation:acetylation contact ratio', 'Normalized acetylation:methylation contact ratio')
% xlabel('radius (nm)')
% ylabel('Normalized contact ratio')
 meth_domain=mean(meth_within,1)./((4/3)*(pi()*((10:10:sdis)/32).^3));
 ac_domain=mean(ac_within_ac,1)./((4/3)*(pi()*((10:10:sdis)/32).^3));
% figure; plot(10:10:sdis,mean(meth_within,1)./((4/3)*(pi()*((10:10:sdis)/32).^3)),'g')
% hold on; plot(10:10:sdis,mean(ac_within_ac,1)./((4/3)*(pi()*((10:10:sdis)/32).^3)),'r')
% legend('H3K27me3 RDF', 'H3K27ac RDF')
% xlabel('radius (nm)')
% ylabel('RDF')
% xlim([0,300])

meth_RDF=mean(meth_within,1)./((4/3)*(pi()*((10:10:sdis)/32).^3));
ac_RDF=mean(ac_within_ac,1)./((4/3)*(pi()*((10:10:sdis)/32).^3));
meth_ac_RDF=mean(ac_within,1)./((4/3)*(pi()*((10:10:sdis)/32).^3));
ac_meth_RDF=mean(meth_within_ac,1)./((4/3)*(pi()*((10:10:sdis)/32).^3));
% figure; plot(10:10:sdis,(mean(meth_within,1)./((4/3)*(pi()*((10:10:sdis)/32).^3)))./min((mean(meth_within,1)./((4/3)*(pi()*((10:10:sdis)/32).^3)))),'g')
% hold on; plot(10:10:sdis,(mean(ac_within_ac,1)./((4/3)*(pi()*((10:10:sdis)/32).^3)))./min((mean(ac_within_ac,1)./((4/3)*(pi()*((10:10:sdis)/32).^3)))),'r')
% xlim([0,300])
% legend('Normalized H3K27me3 RDF', 'Normalized H3K27ac RDF')
% xlabel('radius (nm)')
% ylabel('Normalized RDF')


% figure; plot(10:10:sdis,meth_RDF./meth_ac_RDF,'g')
% hold on; plot(10:10:sdis,ac_RDF./ac_meth_RDF,'r')
%%
%make nanodomains then do contact ratio
px=10;
xm=meth(:,3); ym=meth(:,4); zm=meth(:,12);
xa=acet(:,3); ya=acet(:,4); za=acet(:,12);

minx=min([xm;xa]); miny=min([ym;ya]); minz=min([zm;za]);
xm=(xm+.01-minx)/px; ym=(ym+.01-miny)/px; zm=(zm+.01-(minz))/px;%(pixel size will be 8nm; you can choose other number)
xa=(xa+.01-minx)/px; ya=(ya+.01-miny)/px; za=(za+.01-(minz))/px;%(pixel size will be 8nm; you can choose other number)
xcomb=[xm;xa]; ycomb=[ym;ya]; zcomb=[zm;za];
imeth=zeros(floor(max(xcomb)+5),floor(max(ycomb)+5),floor(max(zcomb)+5));
meth_hmap=imeth;
imeth_thresh=repmat((imeth(:,:,1)),[1,1,3]);
a1=length(xm);
[meth_dens_sort,order_m]=sort(pix_meth_in_meth);
[ac_dens_sort,order_a]=sort(pix_ac_in_ac);
low_cutoff_meth=round(.05*length(meth_dens_sort));
meth_dens_sort((end-low_cutoff_meth):end)=0;
low_cutoff_ac=round(.05*length(ac_dens_sort));
ac_dens_sort((end-low_cutoff_ac):end)=0;
ac_mask=(ac_dens_sort)>0;
meth_mask=(meth_dens_sort)>0;
xm_thresh=nonzeros(xm(order_m).*meth_mask');
ym_thresh=nonzeros(ym(order_m).*meth_mask');
zm_thresh=nonzeros(zm(order_m).*meth_mask');

xa_thresh=nonzeros(xa(order_a).*ac_mask');
ya_thresh=nonzeros(ya(order_a).*ac_mask');
za_thresh=nonzeros(za(order_a).*ac_mask');
a3=length(xm_thresh);
for ind=1:a1
   imeth(floor(xm(ind))+1,floor(ym(ind))+1,floor(zm(ind))+1)=imeth(floor(xm(ind))+1,floor(ym(ind))+1,floor(zm(ind))+1)+1;
   meth_hmap(floor(xm(ind))+1,floor(ym(ind))+1,floor(zm(ind))+1)=imeth(floor(xm(ind))+1,floor(ym(ind))+1,floor(zm(ind))+1)+(pix_meth_in_ac(ind)./pix_meth_in_meth(ind));
end

%Generate the image of the second split
imac=zeros(floor(max(xcomb)+5),floor(max(ycomb)+5),floor(max(zcomb)+5));
ac_hmap=imac;
imac_thresh=repmat((imac(:,:,1)),[1,1,3]);
a2=length(xa);
a4=length(xa_thresh);
for ind=1:a2
   imac(floor(xa(ind))+1,floor(ya(ind))+1,floor(za(ind))+1)=imac(floor(xa(ind))+1,floor(ya(ind))+1,floor(za(ind))+1)+1;
   ac_hmap(floor(xa(ind))+1,floor(ya(ind))+1,floor(za(ind))+1)=imac(floor(xa(ind))+1,floor(ya(ind))+1,floor(za(ind))+1)+(pix_ac_in_meth(ind)./pix_ac_in_ac(ind)); 

end


%%

%%
cell_img2bin = imbinarize(imgaussfilt(sum(imeth+imac,3),1), 0.17);
cell_img2bin = imgaussfilt(double(cell_img2bin),2);
cell_img2bin=imfill(cell_img2bin,'holes');
cell_img2bin2 = bwareaopen(logical(cell_img2bin),(500000));
cell_img2bin=cell_img2bin2.*cell_img2bin;


%cell_img2bin= ~bwareaopen(~cell_img2bin, 80000);
figure; imagesc(logical(cell_img2bin));
objects=bwconncomp(cell_img2bin);
imeth2=sum(imeth,3);
imac2=sum(imac,3);
for n=1
    blank=zeros(size(cell_img2bin,1),size(cell_img2bin,2));
    blank(objects.PixelIdxList{1,n})=1;
    blank2=zeros(size(cell_img2bin,1),size(cell_img2bin,2),3);
    se = strel('disk',5);
    new_im=blank;
    index=0;
    recon=0;
    while sum(sum(new_im))>0
        if index==95
            fddd=1;
        end
        new_imt=new_im;
        new_im=imerode(new_im,se);
        dif=new_imt-new_im;
        sarea=sum(sum(dif))*(px/1000)*(px/1000);
        methyl=sum(sum(dif.*imeth2));
        acyl=sum(sum(dif.*imac2));
        if sum(sum(dif.*imeth2))>85
            recon=1;
   
        else
            recon=0;
        end
        blank2(:,:,1)=dif.*imac2;
        blank2(:,:,2)=dif.*imeth2;
        if recon==1
            index=index+1;
            ratio.n(index)=methyl/acyl;
            mdens(index)=methyl/sarea;
            adens(index)=acyl/sarea;
            gifblank1(:,:,index)=dif;
           % gifmeth1(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,2),2)',2)),1);
           % gifac1(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,1),2)',2)),1);
        end
    end
end
% 
% for n=2
%     blank=zeros(size(cell_img2bin,1),size(cell_img2bin,2));
%     blank(objects.PixelIdxList{1,n})=1;
%     blank2=zeros(size(cell_img2bin,1),size(cell_img2bin,2),3);
%     se = strel('disk',5);
%     new_im=blank;
%     index=0;
%     recon=0;
%     while sum(sum(new_im))>0
%         
%         new_imt=new_im;
%         new_im=imerode(new_im,se);
%         dif=new_imt-new_im;
%         sarea=sum(sum(dif))*(px/1000)*(px/1000);
%         methyl=sum(sum(dif.*imeth2));
%         acyl=sum(sum(dif.*imac2));
%         if sum(sum(dif.*imeth2))>85
%             recon=1;
%         end
%         blank2(:,:,1)=dif.*imac2;
%         blank2(:,:,2)=dif.*imeth2;
%         if recon==1
%             index=index+1;
%             ratio.m(index)=methyl/acyl;
%             gifblank2(:,:,index)=dif;
%             gifmeth2(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,2),2)',2)),1);
%             gifac2(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,1),2)',2)),1);
%         end
%     end
% end
% 
% %%
 filename='ring';
 filename4=[filename,'_binary1','.tif']; % Specify the output file name
 [~,~,d]=size(gifblank1);
 for idx = 1:d
     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
     if idx == 1
         imwrite(uint8(gifblank1(:,:,idx)),filename4);
     else
         imwrite(uint8(gifblank1(:,:,idx)),filename4,'WriteMode','append');
     end
 end
% 
% filename3=[filename,'_meth1','.tif']; % Specify the output file name
% [~,~,d]=size(gifmeth1);
% gifmeth1=gifmeth1*255;
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(uint8(gifmeth1(:,:,idx)),filename3);
%     else
%         imwrite(uint8(gifmeth1(:,:,idx)),filename3,'WriteMode','append');
%     end
% end
% 
% filename2=[filename,'_ac1','.tif']; % Specify the output file name
% [~,~,d]=size(gifac1);
% gifac1=gifac1*255;
% 
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(uint8(gifac1(:,:,idx)),filename2);
%     else
%         imwrite(uint8(gifac1(:,:,idx)),filename2,'WriteMode','append');
%     end
% end

% filename1=[filename,'_binary2','.tif']; % Specify the output file name
% [~,~,d]=size(gifblank2);
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(uint8(gifblank2(:,:,idx)),filename1);
%     else
%         imwrite(uint8(gifblank2(:,:,idx)),filename1,'WriteMode','append');
%     end
% end
% 
% filename5=[filename,'_meth2','.tif']; % Specify the output file name
% [~,~,d]=size(gifmeth2);
% gifmeth2=gifmeth2*255;
% 
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(uint8(gifmeth2(:,:,idx)),filename5);
%     else
%         imwrite(uint8(gifmeth2(:,:,idx)),filename5,'WriteMode','append');
%     end
% end
% 
% filename6=[filename,'_ac2','.tif']; % Specify the output file name
% [~,~,d]=size(gifac2);
% gifac2=gifac2*255;
% 
% for idx = 1:d
%     %[A,map] = rgb2ind(temp2(:,:,:,idx),256);
%     if idx == 1
%         imwrite(uint8(gifac2(:,:,idx)),filename6);
%     else
%         imwrite(uint8(gifac2(:,:,idx)),filename6,'WriteMode','append');
%     end
% end
ratio=ratio.n;
%%
for n1=1:length(xm_thresh)
    if cell_img2bin(floor(xm_thresh(n1))+1,floor(ym_thresh(n1))+1)==0
        xm_thresh(n1)=nan;
        ym_thresh(n1)=nan;
    end
end
for n2=1:length(xa_thresh)
    if cell_img2bin(floor(xa_thresh(n2))+1,floor(ya_thresh(n2))+1)==0
        xa_thresh(n2)=nan;
        ya_thresh(n2)=nan;
    end
end
xm_thresh(isnan(xm_thresh))=[];
ym_thresh(isnan(ym_thresh))=[];
xa_thresh(isnan(xa_thresh))=[];
ya_thresh(isnan(ya_thresh))=[];
%%
sarea=sum(sum(logical(cell_img2bin))).*((px/1000).^2);
meth_dens=(length(xm)/sarea);
ac_dens=(length(xa)/sarea);
%factorm=-0.02372*meth_dens+58.4585;
load('methyl_epsilon_surround2')
load('acetyl_epsilon_surround2')
factorm=coef*(meth_dens.^exp);
factora=coefa*(ac_dens.^expa);
%factora=-0.02372*ac_dens+58.4585;

idxm=dbscan([xm_thresh,ym_thresh],factorm/px,8);
idxa=dbscan([xa_thresh,ya_thresh],factora/px,8);
 figure;
for nm=1:max(idxm)
    idxm_t=find(idxm==nm);
    xm_t=xm_thresh(idxm_t);
    ym_t=ym_thresh(idxm_t);
    zm_t=zm_thresh(idxm_t);
     hold on; plot(xm_t,ym_t,'.')
    meth_cent(nm,:)=[mean(xm_t),mean(ym_t),mean(zm_t)];
    color=rand(1,1,3);
    color=color/sum(color);
    for ind=1:length(xm_t)
        imeth_thresh(floor(xm_t(ind))+1,floor(ym_t(ind))+1,:)=color;
    end
end
 figure;
for na=1:max(idxa)
    idxa_t=find(idxa==na);
    xa_t=xa_thresh(idxa_t);
    ya_t=ya_thresh(idxa_t);
    za_t=za_thresh(idxa_t);
     hold on; plot(xa_t,ya_t,'.')
     ac_cent(na,:)=[mean(xa_t),mean(ya_t),mean(za_t)];
     color=rand(1,1,3);
     color=color/sum(color);
     for ind=1:length(xa_t)
        
         imac_thresh(floor(xa_t(ind))+1,floor(ya_t(ind))+1,:)=color;
    end
end


%meth_peaks=
%%
%Find contact ratio of nanodomains
xm=xm*px;
xa=xa*px;

ym=ym*px;
ya=ya*px;
zm=zm*px;
za=za*px;
meth_cent=meth_cent*px;
ac_cent=ac_cent*px;

[xm_c,order]=sort(meth_cent(:,1)); ym_c=meth_cent(order,2); zm_c=meth_cent(order,3);
[xa_c,order]=sort(ac_cent(:,1)); ya_c=ac_cent(order,2); za_c=ac_cent(order,3);


meth_withinc=zeros(length(xm_c),25);
ac_withinc=zeros(length(xm_c),25);
sdis=1000;
for ind2=1:length(xm_c)
       mask=(abs(xm-xm_c(ind2))<=sdis).*(abs(ym-ym_c(ind2))<=sdis);
       mask=logical(mask); 
      
       %%
       [ind,dist]=rangesearch([xm(mask),ym(mask)],[xm_c(ind2),ym_c(ind2)],sdis);
       ztemp=zm(mask);
       ztemp=ztemp(ind{1});
       distz=abs(ztemp-zm_c(ind2));
       
       for nnn=10:10:sdis
            meth_withinc(ind2,nnn/10)=sum((dist{1}<=nnn))-1;
            tdistz=abs(distz(1:length(nonzeros((dist{1}<=(nnn-10))))));
            mdistz(ind2,nnn/10)=nanmean(abs((distz((length(tdistz)+1):(length(nonzeros((dist{1}<=(nnn)))))))));
            
        end
        
        mask=(abs(xa-xm_c(ind2))<=sdis).*(abs(ya-ym_c(ind2))<=sdis).*(abs(za-zm_c(ind2))<=sdis);
        mask=logical(mask);    
       [ind,dist]=rangesearch([xa(mask),ya(mask),za(mask)],[xm_c(ind2),ym_c(ind2),zm_c(ind2)],sdis);
        for nnn=10:10:sdis
            ac_withinc(ind2,nnn/10)=sum((dist{1}<=nnn));
            
        end
end

meth_within_ac_c=zeros(length(xa_c),25);
ac_within_ac_c=zeros(length(xa_c),25);
for ind2=1:length(xa_c)
       mask=(abs(xm-xa_c(ind2))<=sdis).*(abs(ym-ya_c(ind2))<=sdis);
       mask=logical(mask);
       [ind,dist]=rangesearch([xm(mask),ym(mask),zm(mask)],[xa_c(ind2),ya_c(ind2),za_c(ind2)],sdis);
      
       for nnn=10:10:sdis
            meth_within_ac_c(ind2,nnn/10)=sum((dist{1}<=nnn));

        end
        
        mask=(abs(xa-xa_c(ind2))<=sdis).*(abs(ya-ya_c(ind2))<=sdis);
        mask=logical(mask);
       [ind,dist]=rangesearch([xa(mask),ya(mask)],[xa_c(ind2),ya_c(ind2)],sdis);
       %[indz,distz]=rangesearch([za(mask)],[za_c(ind2)],sdis);
         ztemp=za(mask);
       ztemp=ztemp(ind{1});
       distz=abs(ztemp-za_c(ind2));
       
       for nnn=10:10:sdis
            ac_within_ac_c(ind2,nnn/10)=sum((dist{1}<=nnn))-1;
            tdistz=abs(distz(1:length(nonzeros((dist{1}<=(nnn-10))))));
            adistz(ind2,nnn/10)=nanmean(abs((distz((length(tdistz)+1):(length(nonzeros((dist{1}<=(nnn)))))))));
       end
end
%%
%%
% figure; plot((10:10:sdis),mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1),'g')
% hold on; plot((10:10:sdis),mean(ac_within_ac_c(:,:),1)./mean(meth_within_ac_c(:,:),1),'r')
% legend('methylation:acetylation contact ratio', 'acetylation:methylation contact ratio')
% xlabel('radius (nm)')
% ylabel('contact ratio')


% figure; plot((10:10:sdis),(mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1))./(max((mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1)))),'g')
% hold on; plot((10:10:sdis),(mean(ac_within_ac_c(:,:),1)./mean(meth_within_ac_c(:,:),1))./(max((mean(ac_within_ac_c(:,:),1)./mean(meth_within_ac_c(:,:),1)))),'r')
meth_ac_con_c=mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1);
ac_meth_con_c=mean(ac_within_ac_c(:,:),1)./mean(meth_within_ac_c(:,:),1);
% figure; plot((10:10:sdis),(mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1))./meth_ac_con_c(end),'g')
% hold on; plot((10:10:sdis),(mean(ac_within_ac_c(:,:),1)./mean(meth_within_ac_c(:,:),1))./ac_meth_con_c(end),'r')
% legend('Methylation:acetylation nanodomain centroid contact ratio','Acetylation:Methylation nanodomain contact ratio')
% xlabel('radius (nm)')
% ylabel('Normalized contact ratio')
%%
%save('well4_ALDH+_polyco_ch2_img010','sdis','meth_within','ac_within','ac_within_ac','meth_within_ac','meth_ac_con','ac_meth_con','meth_domain','ac_domain','meth_ac_con_c','ac_meth_con_c','meth_withinc','ac_withinc','ac_within_ac_c','meth_within_ac_c')
%%
for ind=1:max(idxm)
indices=find(idxm==ind);
x_diam=max(xm_thresh(indices))-min(xm_thresh(indices));
y_diam=max(ym_thresh(indices))-min(ym_thresh(indices));
diam_m(ind)=(x_diam+y_diam)/2;
end

for ind=1:max(idxa)
indices=find(idxa==ind);
x_diam=max(xa_thresh(indices))-min(xa_thresh(indices));
y_diam=max(ya_thresh(indices))-min(ya_thresh(indices));
diam_a(ind)=(x_diam+y_diam)/2;
end
%%
meth_ab=length(xm);
ac_ab=length(xa);
%%

save(simnames{k},'cell_img2bin','sdis','C1_within','C2_within','C2_within_C2','C1_within_C2','C1_C2_con','C2_C1_con','C1_domain','C2_domain','C1_C2_con_c','C2_C1_con_c','C1_withinc','C2_withinc','C2_within_C2_c','C1_within_C2_c','ratio','diam_m','diam_a','mdens','adens','adistz','mdistz')
%%
clear;
end
