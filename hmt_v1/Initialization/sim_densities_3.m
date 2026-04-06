
clear
load('sim_param')
nNucleoli=2;
px=10;

nucleus=sim_nucleus_3d(nNucleoli);
nucleus=nucleus(:,:,1:length(N));
figure; imagesc(nucleus(:,:,1))
name={'simdens1';'simdens2';'simdens3';'simdens4';'simdens5';'simdens6';'simdens7';'simdens8';'simdens9';'simdens10'};
%%

for n=1:10
    
    sim_histones_3d_2([name{n},'surround2'],nmeth,nac,ratiom*(n),ratioa*(n),ratiomz*(n),ratioaz*(n),nseedm,nseeda,sa,N,frame,meth_between,ac_between_ac,nucleus)
end
% 
% for n=1:10
%     
%     sim_histones_3d_2([name{n},'seed'],nmeth,nac,ratcool,ratcoola,ratcoolz,ratcoolaz,nseedm*n,nseeda*n,sa,N,frame,meth_between,ac_between_ac,nucleus)
% end
