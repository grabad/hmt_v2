%% match %first order localizations to zeroth order localizations 
%and corresponding ground truth positions
%addpath('\\Foilserver_2\smlm_eb53\Ben\DWP_Tests')
%step 1: run denstest on the biggest dataset
%Step 2: Run Simulation_optimization to find ground truth nanodomain size
%Step 3: Run parameter_optimization to find the ideal epsilon value to
%achieve ground truth nanodomain size
fnam1='well_2_2_me3';
fnam0='well_2_2_ac';
fname0=sprintf('%s.csv',fnam0);
fname1=sprintf('%s.csv',fnam1);
pxsz=110;
% fname1='still_first.csv';
% fname0='still.csv';
TR0=csvread(fname0,1,0);
TR1=csvread(fname1,1,0);
simnames={'dens1','dens2','dens3','dens4','dens5','dens6','dens7','dens8','dens9','dens10'};

for k=1:length(simnames)
TRme3=TR1(1:(11-k):end,:);
TRac=TR0(1:(11-k):end,:);
A6={'id','frame','x [nm]','y [nm]','sigmax [nm]','sigmay [nm]','intensity [photon]','offset [photon]','uncertainty [nm]', 'centroid [nm]','bgstd', 'z [nm]', 'z uncertainty','ydif','centroid2','fwhm0','fwhm1'};
savefile=[simnames{k},'_me3.csv'];
writecell(A6,savefile)
dlmwrite(savefile,TRme3,'delimiter',',','-append');
savefile2=[simnames{k},'_ac.csv'];
writecell(A6,savefile2)
dlmwrite(savefile2,TRac,'delimiter',',','-append');
end
