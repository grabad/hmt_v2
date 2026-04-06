function [a_width,m_width,dens,densa]=find_gt_width(meth_between,nmeth,ac_between_ac,nac)
sdis=200;
for nnn=10:10:sdis
    
[N,edges]= histcounts((meth_between(:,nnn/10)/nmeth)*nmeth,'Normalization','probability');
for n2=1:length(N)
   edgemid=(edges(n2)+edges(n2+1))/2;
   weight(n2)=N(n2)*edgemid;
end
mcount(nnn/10)=sum(weight);
saream(nnn/10)=pi()*(nnn^2)-pi*((nnn-10)^2);
end
dens=mcount./saream;

for nnn=10:10:sdis
    
[Na,edgesa]= histcounts((ac_between_ac(:,nnn/10)/nac)*nac,'Normalization','probability');
for n2=1:length(Na)
   edgemida=(edgesa(n2)+edgesa(n2+1))/2;
   weighta(n2)=Na(n2)*edgemida;
end
mcounta(nnn/10)=sum(weighta);
sareaa(nnn/10)=pi()*(nnn^2)-pi*((nnn-10)^2);
end
densa=mcounta./sareaa;
dens=dens-min(dens);
densa=densa-min(densa);
m_width=find(dens<(0.2*max(dens)));
m_width=m_width(1);m_width=2*((m_width*10));
a_width=find(densa<(0.2*max(densa)));
a_width=a_width(1);a_width=2*((a_width*10));
end