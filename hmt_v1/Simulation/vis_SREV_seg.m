fname={'config-14.txt','config-3.txt','config-4.txt','config-5.txt','config-6.txt','config-7.txt','config-8.txt','config-9.txt','config-10.txt','config-11.txt','config-12.txt','config-13.txt','config-13.txt','config-15.txt','config-16.txt','config-17.txt','config-18.txt','config-19.txt','config-20.txt'};
for file=1:1
fname={'config-14.txt','config-3.txt','config-4.txt','config-5.txt','config-6.txt','config-7.txt','config-8.txt','config-9.txt','config-10.txt','config-11.txt','config-12.txt','config-13.txt','config-13.txt','config-15.txt','config-16.txt','config-17.txt','config-18.txt','config-19.txt','config-20.txt'};
fname2={'config_14_seg_alph1p1.xyz','config_3_seg_alph1p1.xyz','config_4_seg_alph1p1.xyz','config_5_seg_alph1p1.xyz','config_6_seg_alph1p1.xyz','config_7_seg_alph1p1.xyz','config_8_seg_alph1p1.xyz','config_9_seg_alph1p1.xyz','config_10_seg_alph1p1.xyz','config_11_seg_alph1p1.xyz','config_12_seg_alph1p1.xyz','config_13_seg_alph1p1.xyz','config_13_seg_alph1p1.xyz','config_15_seg_alph1p1.xyz','config_16_seg_alph1p1.xyz','config_17_seg_alph1p1.xyz','config_18_seg_alph1p1.xyz','config_19_seg_alph1p1.xyz','config_20_seg_alph1p1.xyz'};
fname3={'config_14_seg_alph1p1.mat','config_3_seg_alph1p1.mat','config_4_seg_alph1p1.mat','config_5_seg_alph1p1.mat','config_6_seg_alph1p1.mat','config_7_seg_alph1p1.mat','config_8_seg_alph1p1.mat','config_9_seg_alph1p1.mat','config_10_seg_alph1p1.mat','config_11_seg_alph1p1.mat','config_12_seg_alph1p1.mat','config_13_seg_alph1p1.mat','config_13_seg_alph1p1.mat','config_15_seg_alph1p1.mat','config_16_seg_alph1p1.mat','config_17_seg_alph1p1.mat','config_18_seg_alph1p1.mat','config_19_seg_alph1p1.mat','config_20_seg_alph1p1.mat'};
integration=0; %from 0 to 1
coord=txt2mat(fname{file});
xpos=coord(:,2);
ypos=coord(:,3);
zpos=coord(:,4);
test=zeros(length(xpos),1);
mdist=zeros(length(xpos),1);
tic
for n=1:length(xpos)
rngx=max(xpos)-min(xpos);
rngy=max(ypos)-min(ypos);
rngz=max(zpos)-min(zpos);
fact=70;
maskx=abs((xpos-xpos(n)))<=rngx/fact;
masky=abs((ypos-ypos(n)))<=rngy/fact;
maskz=abs((zpos-zpos(n)))<=rngz/fact;
maskfin=maskx.*masky.*maskz;
while sum(maskfin)<10
    fact=fact/2;
    maskx=abs((xpos-xpos(n)))<=rngx/fact;
    masky=abs((ypos-ypos(n)))<=rngy/fact;
    maskz=abs((zpos-zpos(n)))<=rngz/fact;
    maskfin=maskx.*masky.*maskz;
end
test(n)=sum(maskfin);
[idx,D]=knnsearch([nonzeros((xpos+1000).*maskfin)-1000,nonzeros((ypos+1000).*maskfin)-1000,nonzeros((zpos+1000).*maskfin)-1000],[xpos(n),ypos(n),zpos(n)],'K',11);
mdistn(n)=mean(D(2:end));
if isnan(mdistn(n))
    g=3;
end
end
toc
%%
cutoff=1.165;
mask1=((mdistn(1:end))>cutoff);
masknoise=(rand(length(mask1),1)<(.5*integration));
mask1=abs(mask1-masknoise');
eux=xpos(find(mask1));
euy=ypos(find(mask1));
euz=zpos(find(mask1));
mask2=1-mask1;
hetx=xpos(find(mask2));
hety=ypos(find(mask2));
hetz=zpos(find(mask2));
outputFileName = fname2{file};
%%
% zmaskhet=hetz<0;
% zmaskeu=euz<0;
% eux=eux(zmaskeu);
% euy=euy(zmaskeu);
% euz=euz(zmaskeu);
% 
% hetx=hetx(zmaskhet);
% hety=hety(zmaskhet);
% hetz=hetz(zmaskhet);

%%
gen_ovito_fin(eux,euy,euz,hetx,hety,hetz,outputFileName)
save(fname3{file},'eux','euy','euz','hetx','hety','hetz','xpos','ypos','zpos','mdistn')
%%
% figure; plot3(eux,euy,euz,'.r','MarkerSize',5)
% hold on; plot3(hetx,hety,hetz,'.g','MarkerSize',5)
%%
clear;
end