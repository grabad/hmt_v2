load('well_2_1.mat')
[surroundme1,surroundac1,meth_pos1,ac_pos1,meth_n1,ac_n1,meth_ab1,ratio1,diam_m1,diam_a1,adens1,mdens1,no_meth1,no_ac1,nclusm1,nclusa1]=compile_function(ac_meth_con,ac_within_ac,diam_a,diam_m,meth_ac_con,meth_within,ratio,sdis,mdens,adens);
sa1=sum(sum(logical(cell_img2bin))).*(.01^2); mpeak1=max(meth_ac_con)/mean(meth_ac_con(80:100));apeak1=max(ac_meth_con)/mean(ac_meth_con(80:100));
save('well_2','meth_pos1','ac_pos1','meth_n1','ac_n1','meth_ab1','ratio1','diam_m1','diam_a1','adens1','mdens1','no_meth1','no_ac1','sa1','nclusm1','nclusa1');

load('well_2_2.mat')
[surroundme2,surroundac2,meth_pos2,ac_pos2,meth_n2,ac_n2,meth_ab2,ratio2,diam_m2,diam_a2,adens2,mdens2,no_meth2,no_ac2,nclusm2,nclusa2]=compile_function(ac_meth_con,ac_within_ac,diam_a,diam_m,meth_ac_con,meth_within,ratio,sdis,mdens,adens);
sa2=sum(sum(logical(cell_img2bin))).*(.01^2); mpeak2=max(meth_ac_con)/mean(meth_ac_con(80:100));apeak2=max(ac_meth_con)/mean(ac_meth_con(80:100));
save('well_2','meth_pos2','ac_pos2','meth_n2','ac_n2','meth_ab2','ratio2','diam_m2','diam_a2','adens2','mdens2','no_meth2','no_ac2','sa2','nclusm2','nclusa2','-append');


no_meth=[no_meth1;no_meth2];
no_ac=[no_ac1;no_ac2];
sa=[sa1;sa2];

ratio=[ratio1;ratio2];
ac_n=[ac_n1;ac_n2];
meth_n=[meth_n1;meth_n2];
figure; plot(10:10:200,mean(ac_n,1),'r','LineWidth',4)
hold on; plot(10:10:200,mean(meth_n,1),'g','LineWidth',4)
legend('acetylation','methylation')
ylabel('mean self-nonself contact')
xlabel('distance (nm)')
ylim([0,1])

meth_pos=[meth_pos1;meth_pos2];
ac_pos=[ac_pos1;ac_pos2];
abundance=[meth_ab1;meth_ab2];
mdens=[mdens1;mdens2];
adens=[adens1;adens2];
surroundme=mean([surroundme1;surroundme2]);
surroundac=mean([surroundac1;surroundac2]);
densitym=no_meth./sa;
densitya=no_ac./sa;
mpeak=[mpeak1;mpeak2];
apeak=[apeak1;apeak2];

mean_ratio=mean(ratio,1);
madens=mean(adens,1);
mmdens=mean(mdens,1);
figure; plot(95:-1:5,mean_ratio(5:95)./max(mean_ratio(5:95)),'b','LineWidth',4);
ylim([0.2,1])
ylabel('mean ratio methylation:acetylation')
xlabel('Radius (%)')
figure; plot(95:-1:5,madens(5:95)/max(madens(5:95)),'r','LineWidth',4);
ylim([0.2,1])
hold on; plot(95:-1:5,mmdens(5:95)/max(mmdens(5:95)),'g','LineWidth',4);
ylim([0.2,1])

diam_methylation=10*[diam_m1;diam_m2];
diam_acetylation=10*[diam_a1;diam_a2];

nclusm=[nclusm1;nclusm2];
nclusa=[nclusa1;nclusa2];


save('well_2','ac_n','meth_n','meth_pos','ac_pos','diam_methylation','diam_acetylation','mean_ratio','abundance','mdens','adens','no_meth','no_ac','nclusm','nclusa','-append');
