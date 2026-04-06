function [surroundme1,surroundac1,meth_pos1,ac_pos1,meth_n1,ac_n1,meth_ab1,ratio1,diam_m1,diam_a1,adens1,mdens1,no_meth1,no_ac1,nclusm1,nclusa1]=compile_function(ac_meth_con,ac_within_ac,diam_a,diam_m,meth_ac_con,meth_within,ratio,sdis,mdens,adens)
sdis=10:10:sdis;
surroundme1=mean(meth_within)./length(meth_within);
surroundac1=mean(ac_within_ac)./length(ac_within_ac);
no_meth1=size(meth_within,1);
no_ac1=size(ac_within_ac,1);
nclusa1=length(diam_a);
nclusm1=length(diam_m);
minmeth=min(meth_ac_con); maxmeth=max(meth_ac_con);
minac=min(ac_meth_con); maxac=max(ac_meth_con);
meth_n=(meth_ac_con(1:20)-minmeth)./(maxmeth-minmeth);
ac_n=(ac_meth_con(1:20)-minac)./(maxac-minac);
meth_n1=meth_n;
ac_n1=ac_n;
meth_ab=length(meth_within)./length(ac_within_ac);
sdis=sdis(1:20);
IND=1;
while ac_n(IND+1)>=ac_n(IND)
    IND=IND+1
end
Istart=IND;
while ((IND+1)<20)&&(ac_n(IND+1)<ac_n(IND))
    IND=IND+1;
end
Iend=IND+1;
ac_n=ac_n(Istart:Iend);
sdisa=sdis(Istart:Iend);
IND=1;
while meth_n(IND+1)>=meth_n(IND)
    IND=IND+1
end
Istart=IND;
while ((IND+1)<20)&&(meth_n(IND+1)<meth_n(IND))
    IND=IND+1;
end
Iend=IND;
meth_n=meth_n(Istart:Iend);
sdism=sdis(Istart:Iend);
meth_f = fit(sdism',meth_n','exp1');
ac_f = fit(sdisa',ac_n','exp1');
meth_pos1=meth_f.b;
ac_pos1=ac_f.b;
ratio1=resample(ratio,100,length(ratio));
meth_ab1=meth_ab;
diam_m1=mean(diam_m);
diam_a1=mean(diam_a);
mdens1=resample(mdens,100,length(mdens));
adens1=resample(adens,100,length(adens));

end