function [ratiom_fin,ratioa_fin,ratiomz_fin,ratioaz_fin,meth_between,ac_between_ac,dif]=create_variables(meth_withinc,mdistz,ac_within_ac_c,adistz,sim_meth_withinc,sim_mdistz,sim_ac_within_ac_c,sim_adistz,ratiom_i,ratioa_i,ratiomz_i,ratioaz_i)
%%

sim_meth_between(:,1)=sim_meth_withinc(:,1);
for n=size(meth_withinc,2):-1:2
    sim_meth_between(:,n)=sim_meth_withinc(:,n)-sim_meth_withinc(:,n-1);
end
sim_meth_between(sim_meth_between<0)=0;

sim_ac_between_ac(:,1)=sim_ac_within_ac_c(:,1);
for n=size(sim_ac_within_ac_c,2):-1:2
    sim_ac_between_ac(:,n)=sim_ac_within_ac_c(:,n)-sim_ac_within_ac_c(:,n-1);
end
sim_ac_between_ac(sim_ac_between_ac<0)=0;
%%
%%
for n=1:1:20
%figure; histogram(sim_meth_between(:,n));
xxx(n)=mean(sim_meth_between(:,n));
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');
end

for n=1:1:20
%figure; histogram(sim_ac_between_ac(:,n));
xxxa(n)=mean(sim_ac_between_ac(:,n));
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');
end

%figure; histogram(sim_meth_betweenz(:,n));
xxxz=nanmean(sim_mdistz);
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');


%figure; histogram(sim_ac_between_acz(:,n));
xxxaz=nanmean(sim_adistz);
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');
%%

meth_between(:,1)=meth_withinc(:,1);
for n=size(meth_withinc,2):-1:2
    meth_between(:,n)=meth_withinc(:,n)-meth_withinc(:,n-1);
end
meth_between(meth_between<0)=0;

ac_between_ac(:,1)=ac_within_ac_c(:,1);
for n=size(ac_within_ac_c,2):-1:2
    ac_between_ac(:,n)=ac_within_ac_c(:,n)-ac_within_ac_c(:,n-1);
end
ac_between_ac(ac_between_ac<0)=0;
%%


%%
%norm_factm=length(meth_within);
for n=1:1:20
%figure; histogram(meth_between(:,n));
xx(n)=mean(meth_between(:,n));
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');
end

for n=1:1:20
%figure; histogram(ac_between_ac(:,n));
xxa(n)=mean(ac_between_ac(:,n));
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');
end

%figure; histogram(meth_betweenz(:,n));
xxz=nanmean(mdistz);
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');


%figure; histogram(ac_between_acz(:,n));
xxaz=nanmean(adistz);
%[N(n,:),edges(n,:)] = histcounts(meth_within(:,n)/norm_factm,'Normalization','probability');

%%
ratiom=xxx./xx;

 ratioa=xxxa./xxa;

 ratiomz=xxxz./xxz;

 ratioaz=xxxaz./xxaz;
% load('fit_res_fin_3d.mat')
  ratiom_fin=ratiom_i.*ratiom;
  ratioa_fin=ratioa_i.*ratioa;
    ratiomz_fin=ratiomz_i(1:20).*ratiomz(1:20);
  ratioaz_fin=ratioaz_i(1:20).*ratioaz(1:20);
  dif=max(abs(([ratiom,ratioa,ratiomz,ratioaz])-1));
%   ratcool_fin=ratcool;
%   ratcoola_fin=ratcoola;

%save('fit_res_fin_3d','ratcool_fin','ratcoola_fin')
end