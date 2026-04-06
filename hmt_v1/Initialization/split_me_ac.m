function split_me_ac
k=1;
simnames={'Simulation'};
TR0 = csvread([simnames{k},'.csv'],1,0);

maskme=TR0(:,6)<680;
for n=1:size(TR0,2)
    TRme3(:,n)=TR0(maskme,n);
end
maskac=TR0(:,6)>680;
for n=1:size(TR0,2)
    TRac(:,n)=TR0(maskac,n);
end

A6={'id','frame','x [nm]','y [nm]','z [nm]','centroid [nm]'};
writecell(A6,[simnames{k},'_me3.csv'])
dlmwrite([simnames{k},'_me3.csv'],TRme3,'delimiter',',','-append');

A6={'id','frame','x [nm]','y [nm]','z [nm]','centroid [nm]'};
writecell(A6,[simnames{k},'_ac.csv'])
dlmwrite([simnames{k},'_ac.csv'],TRac,'delimiter',',','-append');
clear;

end