function []=sim_histones_3d(name,nmeth,nac,ratcool,ratcoola,ratcoolz,ratcoolaz,nseedm,nseeda,sa,N,frame,meth_between,ac_between_ac)
nNucleoli=2;
px=10;

nucleus=sim_nucleus_3d(nNucleoli);
nucleus=nucleus(:,:,1:length(N));
figure; imagesc(nucleus(:,:,1))
%%
filename = name;
filename2=[filename,'.tif']; % Specify the output file name
nucleus2=uint16(nucleus);
sa_sim=(sum(sum(logical(sum(nucleus,3))))).*(px/1000).^2;
[~,~,d]=size(nucleus);
for idx = 1:d
%[A,map] = rgb2ind(temp2(:,:,:,idx),256);
if idx == 1
imwrite(nucleus2(:,:,idx),filename2);
else
imwrite(nucleus2(:,:,idx),filename2,'WriteMode','append');
end
end
%%
for zsiz=1:size(nucleus,3)
objects=bwconncomp(nucleus(:,:,zsiz));
% imeth2=sum(imeth,3);
% imac2=sum(imac,3);
load('mean_ratio.mat')
load('mean_adens');
load('mean_mdens');

    blank=zeros(size(nucleus,1),size(nucleus,2));
    blank(objects.PixelIdxList{1,1})=1;
    se = strel('disk',5);
    new_im=blank;
    new_im=imfill(new_im);
    index=0;
    recon=0;
    while sum(sum(new_im))>0
        index=index+1;
        new_imt=new_im;
        new_im=imerode(new_im,se);
        dif(:,:,index)=new_imt-new_im;
        
 
            
         
           % gifblank1(:,:,index)=dif;
           % gifmeth1(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,2),2)',2)),1);
           % gifac1(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,1),2)',2)),1);
    end
mean_ratio=resample(mean_ratio,index,length(mean_ratio));
mean_ratio=mean_ratio.^1.5;
mdens=resample(mdens,index,length(mdens));

[a,b,c]=size(dif);
densities_eus=zeros(a,b,c);
for n=1:index
   densities_eus(:,:,n)=dif(:,:,n).*mdens(n);
   dif(:,:,n)=dif(:,:,n).*mean_ratio(n); 
end
dif=sum(dif,3);
densities_eus=sum(densities_eus,3);
densities_eusf(:,:,zsiz)=densities_eus;
diff(:,:,zsiz)=dif;
end


%%
for zsiz2=1:size(nucleus,3)
nucleoli(:,:,zsiz2) = bwareaopen(~nucleus(:,:,zsiz2), 2000);
BW2 = bwareafilt(nucleoli(:,:,zsiz2), 1);
nucleusf= imfill(nucleus(:,:,zsiz2),'holes');
%nucleus = ~bwareaopen(~nucleus, 20);
nucleoli(:,:,zsiz2)=nucleoli(:,:,zsiz2)-BW2;
%figure; imagesc(logical(nucleus(:,:,zsiz2)));
objects=bwconncomp(nucleus(:,:,zsiz2));


    blank=zeros(size(nucleusf,1),size(nucleusf,2));
    mem=blank;
    blank(objects.PixelIdxList{1,1})=1;
    se = strel('disk',5);
    outer=logical(nucleusf)-imerode(logical(nucleusf),se);
    new_im=nucleoli(:,:,zsiz2);
    index=1;
    dif2=new_im;
    while sum(sum((dif2(:,:,index)-mem).*logical(nucleusf)))>0
        if index==95
            fddd=1;
        end
        index=index+1;
        new_imt=new_im;
        new_im=imdilate(new_im,se);
        dif2(:,:,index)=new_im-new_imt;
        %mem=mem+outer.*dif2(:,:,index);

        dif2(:,:,index)=dif2(:,:,index).*logical(nucleusf);
        dif2(:,:,index)=dif2(:,:,index);%+mem;



      

           % gifmeth1(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,2),2)',2)),1);
           % gifac1(:,:,index)=imgaussfilt(double(downsample(downsample(blank2(:,:,1),2)',2)),1);
      
    end
load('ratio')
load('mdens2')

mdens=resample(mdens,index,length(mdens));
ratio=resample(ratio,index,length(ratio));
ratio=ratio/2.5+.45;
endpoint=round(index/4);
[a,b,c]=size(dif2);
densities_oli=zeros(a,b,c);
dif2=double(dif2);
for n=1:endpoint
   densities_oli(:,:,n)=dif2(:,:,n).*mdens(n);
   dif2(:,:,n)=dif2(:,:,n).*ratio(n);
   
end
dif2=sum(dif2,3);
densities_oli=sum(densities_oli,3);

dif2(dif2==1)=0;
dif_fin=max(diff(:,:,zsiz2),dif2);
dif_finf(:,:,zsiz2)=dif_fin;
densities_fin=max(densities_eusf(:,:,zsiz2),densities_oli);
densities_fin=densities_fin.*nucleus(:,:,zsiz2);

densities_finf(:,:,zsiz2)=densities_fin;
end
%%
%load('z_dist.mat');
save('sim_param2')
mind=0;
aind=0;
sarat=sa/sa_sim;
for n34=1:size(nucleus,3)
    nuc_ind=find(nucleus(:,:,n34));
    densities_fin=densities_finf(:,:,n34);
    dif_fin=dif_finf(:,:,n34);
for n=1:(sum(sum(nucleus(:,:,n34))))
    
    seedling=rand;

    meth_thresh=nseedm*(densities_fin(nuc_ind(n)))*1/((561000000*sarat)/N(n34));
    ac_thresh=1-(((nseeda*densities_fin(nuc_ind(n)))/dif_fin(nuc_ind(n)))*1/((425110000*sarat)/N(n34)));
    
    if seedling<meth_thresh
        mind=mind+1;
        [I,J]=ind2sub([size(nucleus,1),size(nucleus,2)],nuc_ind(n));
        methpos(mind,1)=(I-1+rand)*px;
        methpos(mind,2)=(J-1+rand)*px;
        methpos(mind,3)=frame(n34);
    elseif seedling>ac_thresh
        aind=aind+1;
        [I,J]=ind2sub([size(nucleus,1),size(nucleus,2)],nuc_ind(n));
        acpos(aind,1)=(I-1+rand)*px;
        acpos(aind,2)=(J-1+rand)*px;
        acpos(aind,3)=frame(n34);
    else
    end
        
end
end
%%

sdis=200; %nm
norm_factm=nmeth;
norm_facta=nac;

Nstab=length(methpos)
for ind2=1:Nstab
    mindo=ind2;
        for nnn=10:10:sdis
            %surroundme(nnn/10);
            [N,edges] = histcounts((meth_between(:,nnn/10)/norm_factm)*nmeth,'Normalization','probability');
            seedm=rand;
           edges=round(edges/(6*ratcool(nnn/10)));
            cd=0;
            
            clear cdf;
            for n2=1:length(N)
                cd=cd+N(n2);
                cdf(n2)=cd;
            end
            nmtemp=find(seedm<cdf);
            nmtemp=nmtemp(1);
            nmt=round(mean([edges(nmtemp),edges(nmtemp+1)]));
            
            for n3=1:nmt
            mind=mind+1;
            dist=nnn+10*(rand);
            angrand=rand*2*pi();
            xdist=dist*cos(angrand);
            ydist=dist*sin(angrand);
            zdist=(2.8*(sign(rand-.5))*sqrt(((nnn+10*(2*(rand)-1)).^2)/2))./ratcoolz(nnn/10);
            methpos(mind,1)=methpos(mindo,1)+xdist;
            methpos(mind,2)=methpos(mindo,2)+ydist;
            methpos(mind,3)=methpos(mindo,3)+zdist;
            end
        end
            
end
Nastab=length(acpos)
for ind2=1:Nastab
    aindo=ind2;
        for nnn=10:10:sdis
            %surroundme(nnn/10);
            [N,edges] = histcounts((ac_between_ac(:,nnn/10)/norm_facta)*nac,'Normalization','probability');
            seeda=rand;
            edges=round(edges/(6*ratcoola(nnn/10)));

            cd=0;
            clear cdf;
            for n2=1:length(N)
                cd=cd+N(n2);
                cdf(n2)=cd;
            end
            natemp=find(seeda<cdf);
            natemp=natemp(1);
            nat=round(mean([edges(natemp),edges(natemp+1)]));
            
            for n3=1:nat
            aind=aind+1;
            dist=nnn+10*rand;
            angrand=rand*2*pi();
            xdist=dist*cos(angrand);
            ydist=dist*sin(angrand);
            zdist=(2.8*(sign(rand-.5))*sqrt(((nnn+10*(2*(rand)-1)).^2)/2))./ratcoolaz(nnn/10);

            acpos(aind,1)=acpos(aindo,1)+xdist;
            acpos(aind,2)=acpos(aindo,2)+ydist;
            acpos(aind,3)=acpos(aindo,3)+zdist;
            end
        end
            
end
%%
nloc=length(acpos)+length(methpos);
frames=ones(nloc,1);
id=1:nloc;
centroid=[670*ones(length(methpos),1);700*ones(length(acpos),1)];

TRsim=[id',frames,[methpos(:,1);acpos(:,1)],[methpos(:,2);acpos(:,2)],[methpos(:,3);acpos(:,3)],centroid];

A6={'id','frame','x [nm]','y [nm]','z [nm]','centroid [nm]'};
%savefile='Simulation_big.csv';
savefile=[filename,'.csv'];
writecell(A6,savefile)

dlmwrite(savefile,TRsim,'delimiter',',','-append');
end

