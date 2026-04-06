function [meth_pos1,ac_pos1]=compile_function_in(ac_meth_con,meth_ac_con)
sdis=10:10:1000;
if ac_meth_con(1)>10^6
    ac_meth_con=ac_meth_con(2:end); 
end
if meth_ac_con(1)>10^6
    meth_ac_con=meth_ac_con(2:end); 
end


minmeth=min(meth_ac_con); maxmeth=max(meth_ac_con);
minac=min(ac_meth_con); maxac=max(ac_meth_con);
meth_n=(meth_ac_con(1:20)-minmeth)./(maxmeth-minmeth);
ac_n=(ac_meth_con(1:20)-minac)./(maxac-minac);
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


end