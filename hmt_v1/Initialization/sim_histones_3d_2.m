function []=sim_histones_3d_2(name,nmeth,nac,ratcool,ratcoola,ratcoolz,ratcoolaz,nseedm,nseeda,sa,N,frame,meth_between,ac_between_ac,nucleus)
%%
filenamet = name;

load('sim_param2')
clear filename2; clear filename
filename=filenamet;
filename2=[filename,'.tif']; % Specify the output file name

%%
%load('z_dist.mat');
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

