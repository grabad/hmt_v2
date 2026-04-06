function splitmeac(fnam)
fnam1=[fnam,'_fin'];

fname1=sprintf('%s.csv',fnam1);
TR1=csvread(fname1,1,0);
mask1=(TR1(:,10)>655.62).*(TR1(:,10)<688.35);
mask2=(TR1(:,10)>704.715).*(TR1(:,10)<737.445);
ind1=find(mask1);
ind2=find(mask2);
for n1=1:length(ind1)
    TR_me3(n1,:)=TR1(ind1(n1),:);
end
for n2=1:length(ind2)
    TR_ac(n2,:)=TR1(ind2(n2),:);
end
A6={'id','frame','x [nm]','y [nm]','sigmax [nm]','sigmay [nm]','intensity [photon]','offset [photon]','uncertainty [nm]', 'centroid [nm]','bgstd', 'z [nm]', 'z uncertainty','ydif','centroid2','fwhm0','fwhm1'};
savefile=[fnam,'_fin_me3.csv'];
writecell(A6,savefile)

dlmwrite(savefile,TR_me3,'delimiter',',','-append');

savefile2=[fnam,'_fin_ac.csv'];
writecell(A6,savefile2)

dlmwrite(savefile2,TR_ac,'delimiter',',','-append');
end