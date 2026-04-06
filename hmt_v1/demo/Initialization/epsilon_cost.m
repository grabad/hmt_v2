function cost = epsilon_cost(factor,xm,ym,zm,gt)
px=10;
idx=dbscan([xm/px,ym/px],factor/px,8);

%meth_peaks=
%%
%Find contact ratio of nanodomains


diam_m=[];
%%
for ind=1:max(idx)
indices=find(idx==ind);
x_diam=max(xm(indices))-min(xm(indices));
y_diam=max(ym(indices))-min(ym(indices));
diam_m(ind)=(x_diam+y_diam)/2;
end
%%
if isempty(diam_m)
    cost=-100
else
cost=mean(diam_m)-gt;
end
if factor>29.99999
   dd=3; 
end
%densm=meth_ab/(sum(sum(logical(cell_img2bin))).*.01^2);
%save([simnames{k},'_orig'],'diam_a','diam_m','meth_dens','ac_dens')
end