simnames={'noise1_seg_alph1p2';'noise2_seg_alph1p2';'noise3_seg_alph1p2';'noise4_seg_alph1p2';'noise5_seg_alph1p2';'noise6_seg_alph1p2';'noise7_seg_alph1p2';'noise8_seg_alph1p2';'noise9_seg_alph1p2';'noise10_seg_alph1p2';'noise11_seg_alph1p2';'noise12_seg_alph1p2';'noise13_seg_alph1p2';'noise14_seg_alph1p2';'noise15_seg_alph1p2';'noise16_seg_alph1p2';'noise17_seg_alph1p2';'noise18_seg_alph1p2';'noise19_seg_alph1p2';'noise20_seg_alph1p2'};
simnames2={'config_1_seg_alph1p2';'config_2_seg_alph1p2';'config_3_seg_alph1p2';'config_4_seg_alph1p2';'config_5_seg_alph1p2';'config_6_seg_alph1p2';'config_7_seg_alph1p2';'config_8_seg_alph1p2';'config_9_seg_alph1p2';'config_10_seg_alph1p2';'config_11_seg_alph1p2';'config_12_seg_alph1p2';'config_13_seg_alph1p2';'config_14_seg_alph1p2';'config_15_seg_alph1p2';'config_16_seg_alph1p2';'config_17_seg_alph1p2';'config_18_seg_alph1p2';'config_19_seg_alph1p2';'config_20_seg_alph1p2'};
for nth=16:length(simnames)
simnames={'noise1_seg_alph1p2';'noise2_seg_alph1p2';'noise3_seg_alph1p2';'noise4_seg_alph1p2';'noise5_seg_alph1p2';'noise6_seg_alph1p2';'noise7_seg_alph1p2';'noise8_seg_alph1p2';'noise9_seg_alph1p2';'noise10_seg_alph1p2';'noise11_seg_alph1p2';'noise12_seg_alph1p2';'noise13_seg_alph1p2';'noise14_seg_alph1p2';'noise15_seg_alph1p2';'noise16_seg_alph1p2';'noise17_seg_alph1p2';'noise18_seg_alph1p2';'noise19_seg_alph1p2';'noise20_seg_alph1p2'};
simnames2={'config_1_seg_alph1p2';'config_2_seg_alph1p2';'config_3_seg_alph1p2';'config_4_seg_alph1p2';'config_5_seg_alph1p2';'config_6_seg_alph1p2';'config_7_seg_alph1p2';'config_8_seg_alph1p2';'config_9_seg_alph1p2';'config_10_seg_alph1p2';'config_11_seg_alph1p2';'config_12_seg_alph1p2';'config_13_seg_alph1p2';'config_14_seg_alph1p2';'config_15_seg_alph1p2';'config_16_seg_alph1p2';'config_17_seg_alph1p2';'config_18_seg_alph1p2';'config_19_seg_alph1p2';'config_20_seg_alph1p2'};
fnameout = [simnames{nth},'.tif'];

px=100; %nm
%dist=10; %nm
lambda=670; %Emission wavelength in nm
NA=1.4; %Numerica aperture
%nblinks=25;
%ndot=9;
%nit=35;
SH=10; %nm
SD=6; %nm/px
dens=1; %Labelling density
numChains=1; %number of microtubule filaments in the image
FOV=5; %Field of view in microns
nframes=25000; %number of frames
anti_length=0; %Length of antibody in nm; double it if you are using primary and secondary
tub_diam=0; %Circular diameter of tubulin (probably should not change)
numVerts=100000; %1 alpha/beta tubulin epitope per 14.49 nm slong length of microtubule
% k2=0.0135; %Chance of moving from on state to triplet state ms^-1
% k3=0.333e-04; %Chance of moving from triplet state to on state ms^-1
% k4= 0.0011;%Chance of moving from on state to bleached state ms^-1
% frame_rate=40; %ms
% zeroth_seg=500; %Desired photon count of zeroth order
% first_seg=1500; %Desired photon count of first order
Ibgp = 10; %Bacground noise of zeroth order (photons)
Ibg = Ibgp*1.5; %Background noise of first order (photons*100)
RN=1; %readout noise (photons*100)
%dist=15;
%persis=1/35;
    %[a,b,xpos,ypos,zpos,spec_hetr]=sim_grid(px,dist);
 %   [tub,tub_viewable,xpos,ypos,zpos,spec_hetr,identity]=SREV(dens,px,numVerts,FOV);
    load(simnames2{nth});
    xpos=xpos*4.9;
    ypos=ypos*4.9;
    zpos=zpos*4.9;
    coordinates=[xpos,ypos,zpos];
    save('coordinates.txt', 'coordinates', '-ascii');
%     %%
%     box_bounds = [-80 80; -80 80; -80 80]; 
%     cvrt_to_lmp(coordinates,'output.lmp',box_bounds);
    %
%     [xpos2,ind]=unique(xpos);
%     [ind,ind2]=sort(ind);
%     xpos=xpos(ind2);
%     ypos=ypos(ind2);
%     zpos=zpos(ind2);
%     identity=identity(ind2);
%     spec_hetr=spec_hetr(ind2);
    zpos=zpos-mean(zpos);
    %[a,b,c,xtemp,ytemp,ztemp,spec_hetr]=sim_npc(px,lambda,n_npc,anti_length,dens);
    cutoff=1.165;

    identity=2-mask1; % 2 is 647, 1 is 680
    [gtx,gty,gtz,f1x,f1y,f1z,identity2]=sim_anti(xpos,ypos,zpos,identity);
    %[data,gt]=blinking_sim_unrealistic(xpos,ypos,zpos,px,lambda,NA,nframes,nit);
    
    [data,gt]=blinking_sim_unrealistic(f1x,f1y,f1z,px,NA,nframes,FOV,identity2,fnameout);
    n10=0;
    n11=0;
    fnameout647=sprintf('data647.tif');
    fnameout680=sprintf('data680.tif');
    for n9=1:length(f1x)
        if identity2(n9)==2
            n10=n10+1;
            f1x647(n10)=f1x(n9);
            f1y647(n10)=f1y(n9);
            f1z647(n10)=f1z(n9);
            identity2_647(n10)=2;
        else
            n11=n11+1;
            f1x680(n11)=f1x(n9);
            f1y680(n11)=f1y(n9);
            f1z680(n11)=f1z(n9);
            identity2_680(n11)=2;
        end
    end
    [data647,gt647]=blinking_sim_unrealistic(f1x647,f1y647,f1z647,px,NA,nframes,FOV,identity2_647,fnameout647);
    [data680,gt680]=blinking_sim_unrealistic(f1x680,f1y680,f1z680,px,NA,nframes,FOV,identity2_680,fnameout680);

    %[gt]=change_data(gt,nblinks, nit, ndot);
  %  
    make_first_order; %Make first order PSF
    ind647=find(gt(:,6)==671);
    ind680=find(gt(:,6)==705);
    gt647=gt(ind647,:);
    gt680=gt(ind680,:);
    [blinking_first_647]=make_real_fist(sptimg647,SH,gt647,FOV,px,NA,nframes); %Make blinking first order
    [blinking_first_680]=make_real_fist(sptimg680,SH,gt680,FOV,px,NA,nframes); %Make blinking first order
    blinking_first=blinking_first_647+blinking_first_680;
    blinking_first=blinking_first(1:size(data,1),1:size(data,2),:);
    
    fnameout=[simnames{nth},'_z.tif'];
    saveastiff(data, fnameout)
    fnameout=[simnames{nth},'_f.tif'];
    saveastiff(blinking_first, fnameout)

    

    
     add_noise; %Add noise to the image   
    fnameout=[simnames{nth},'_zerothh.tif'];
    saveastiff(data_noise, fnameout)
    fnameout=[simnames{nth},'_firsth.tif'];
    saveastiff(first_noise, fnameout)
    

    

%     s1='tubulin_zeroth_unrealistic.tif'; %name of zeroth order image
%     s2='tubulin_first_unrealistic.tif'; %name of first order image
%     make_picsz(s1,data_noise,s1, size(data_noise,3));  %Download tifs of zeroth order data
%     make_picsz(s2,first_noise,s2, size(first_noise,3)); %Download tifs of first order data
    
    clear;
end
%%
px=100; %nm
simnames={'noise1_seg_alph1p2';'noise2_seg_alph1p2';'noise3_seg_alph1p2';'noise4_seg_alph1p2';'noise5_seg_alph1p2';'noise6_seg_alph1p2';'noise7_seg_alph1p2';'noise8_seg_alph1p2';'noise9_seg_alph1p2';'noise10_seg_alph1p2';'noise11_seg_alph1p2';'noise12_seg_alph1p2';'noise13_seg_alph1p2';'noise14_seg_alph1p2';'noise15_seg_alph1p2';'noise16_seg_alph1p2';'noise17_seg_alph1p2';'noise18_seg_alph1p2';'noise19_seg_alph1p2';'noise20_seg_alph1p2'};
for nth=1:length(simnames)
    QD_analyze_A(simnames{nth},px); %Match zeroth order to first order and ground truth
end
%
simnames={'noise1_seg_alph1p2';'noise2_seg_alph1p2';'noise3_seg_alph1p2';'noise4_seg_alph1p2';'noise5_seg_alph1p2';'noise6_seg_alph1p2';'noise7_seg_alph1p2';'noise8_seg_alph1p2';'noise9_seg_alph1p2';'noise10_seg_alph1p2';'noise11_seg_alph1p2';'noise12_seg_alph1p2';'noise13_seg_alph1p2';'noise14_seg_alph1p2';'noise15_seg_alph1p2';'noise16_seg_alph1p2';'noise17_seg_alph1p2';'noise18_seg_alph1p2';'noise19_seg_alph1p2';'noise20_seg_alph1p2'};
for nth=1:length(simnames)
splitmeac(simnames{nth});
end

    
    for n=1:max(gt(:,2))
        hm(n)=sum(gt(:,2)==n);
    end
hope=find(hm==2);
[frame,order]=sort(gt(:,2));
gt=gt(order,:);
for m=1:length(hope)
    coord=find(gt(:,2)==hope(m));
    x1=gt(coord(1),3);
    y1=gt(coord(1),4);
    x2=gt(coord(2),3);
    y2=gt(coord(2),4);
    int1=gt(coord(1),9);
    int2=gt(coord(2),9);
    dist(m)=sqrt((x1-x2)^2+(y1-y2)^2);
    intd(m)=abs(int1-int2);
end
coord=dist>500;
