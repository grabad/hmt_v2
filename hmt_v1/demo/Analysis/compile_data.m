load('k27_k27_thaw000.mat')
[surroundme1,surroundac1,meth_pos1,ac_pos1,meth_n1,ac_n1,meth_ab1,ratio1,diam_m1,diam_a1,adens1,mdens1,no_meth1,no_ac1,nclusm1,nclusa1]=compile_function(ac_meth_con,ac_within_ac,diam_a,diam_m,meth_ac_con,meth_within,ratio,sdis,mdens,adens);
sa1=sum(sum(logical(cell_img2bin))).*(.01^2); mpeak1=max(meth_ac_con)/mean(meth_ac_con(80:100));apeak1=max(ac_meth_con)/mean(ac_meth_con(80:100));
save('k27_k27_thaw','meth_pos1','ac_pos1','meth_n1','ac_n1','meth_ab1','ratio1','diam_m1','diam_a1','adens1','mdens1','no_meth1','no_ac1','sa1','nclusm1','nclusa1');



no_meth=[no_meth1];
no_ac=[no_ac1];
sa=[sa1];

ratio=[ratio1];
ac_n=[ac_n1];
meth_n=[meth_n1];
figure; plot(10:10:200,mean(ac_n,1),'r','LineWidth',4)
hold on; plot(10:10:200,mean(meth_n,1),'g','LineWidth',4)
legend('acetylation','methylation')
ylabel('mean self-nonself contact')
xlabel('distance (nm)')
ylim([0,1])

meth_pos=[meth_pos1];
ac_pos=[ac_pos1];
abundance=[meth_ab1];
mdens=[mdens1];
adens=[adens1];
surroundme=mean([surroundme1]);
surroundac=mean([surroundac1]);
densitym=no_meth./sa;
densitya=no_ac./sa;
mpeak=[mpeak1];
apeak=[apeak1];

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

diam_methylation=10*[diam_m1];
diam_acetylation=10*[diam_a1];

nclusm=[nclusm1];
nclusa=[nclusa1];


save('k27_k27_thaw','ac_n','meth_n','meth_pos','ac_pos','diam_methylation','diam_acetylation','mean_ratio','abundance','mdens','adens','no_meth','no_ac','nclusm','nclusa','-append');
