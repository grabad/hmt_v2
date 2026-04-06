function [meth_withinc,mdistz,ac_within_ac_c,adistz]=find_within(xm,ym,zm,xa,ya,za,px)
idxm=dbscan([xm/px,ym/px],30/px,8);
idxa=dbscan([xa/px,ya/px],30/px,8);
 figure;
for nm=1:max(idxm)
    idxm_t=find(idxm==nm);
    xm_t=xm(idxm_t);
    ym_t=ym(idxm_t);
    zm_t=zm(idxm_t);
     hold on; plot(xm_t,ym_t,'.')
    meth_cent(nm,:)=[mean(xm_t),mean(ym_t),mean(zm_t)];
    color=rand(1,1,3);
    color=color/sum(color);

end
 figure;
for na=1:max(idxa)
    idxa_t=find(idxa==na);
    xa_t=xa(idxa_t);
    ya_t=ya(idxa_t);
    za_t=za(idxa_t);
     hold on; plot(xa_t,ya_t,'.')
     ac_cent(na,:)=[mean(xa_t),mean(ya_t),mean(za_t)];
     color=rand(1,1,3);
     color=color/sum(color);

end


%meth_peaks=
%%
%Find contact ratio of nanodomains


[xm_c,order]=sort(meth_cent(:,1)); ym_c=meth_cent(order,2); zm_c=meth_cent(order,3);
[xa_c,order]=sort(ac_cent(:,1)); ya_c=ac_cent(order,2); za_c=ac_cent(order,3);


meth_withinc=zeros(length(xm_c),25);
ac_withinc=zeros(length(xm_c),25);
meth_withincz=zeros(length(xm_c),25);

sdis=1000;
for ind2=1:length(xm_c)
       mask=(abs(xm-xm_c(ind2))<=sdis).*(abs(ym-ym_c(ind2))<=sdis);
       mask=logical(mask); 
      
       %%
       [ind,dist]=rangesearch([xm(mask),ym(mask)],[xm_c(ind2),ym_c(ind2)],sdis);
       ztemp=zm(mask);
       ztemp=ztemp(ind{1});
       distz=abs(ztemp-zm_c(ind2));
       
       for nnn=10:10:sdis
            meth_withinc(ind2,nnn/10)=sum((dist{1}<=nnn));
            tdistz=abs(distz(1:length(nonzeros((dist{1}<=(nnn-10))))));
            mdistz(ind2,nnn/10)=nanmean(abs((distz((length(tdistz)+1):(length(nonzeros((dist{1}<=(nnn)))))))));
            
        end
        
        mask=(abs(xa-xm_c(ind2))<=sdis).*(abs(ya-ym_c(ind2))<=sdis);
        mask=logical(mask);    
       [ind,dist]=rangesearch([xa(mask),ya(mask)],[xm_c(ind2),ym_c(ind2)],sdis);
        for nnn=10:10:sdis
            ac_withinc(ind2,nnn/10)=sum((dist{1}<=nnn));
            
        end
end
%%
meth_withincz(ind2,nnn/10)=sum((distz<=nnn));
for ind2=1:length(zm_c)
       mask=(abs(xm-xm_c(ind2))<=sdis).*(abs(ym-ym_c(ind2))<=sdis);
       mask=logical(mask); 
      
       
  

end
%%
meth_within_ac_c=zeros(length(xa_c),25);
ac_within_ac_c=zeros(length(xa_c),25);
ac_within_ac_cz=zeros(length(xa_c),25);
for ind2=1:length(xa_c)
  
        
        mask=(abs(xa-xa_c(ind2))<=sdis).*(abs(ya-ya_c(ind2))<=sdis);
        mask=logical(mask);
       [ind,dist]=rangesearch([xa(mask),ya(mask)],[xa_c(ind2),ya_c(ind2)],sdis);
       %[indz,distz]=rangesearch([za(mask)],[za_c(ind2)],sdis);
         ztemp=za(mask);
       ztemp=ztemp(ind{1});
       distz=abs(ztemp-za_c(ind2));
       
       for nnn=10:10:sdis
            ac_within_ac_c(ind2,nnn/10)=sum((dist{1}<=nnn));
            tdistz=abs(distz(1:length(nonzeros((dist{1}<=(nnn-10))))));
            adistz(ind2,nnn/10)=nanmean(abs((distz((length(tdistz)+1):(length(nonzeros((dist{1}<=(nnn)))))))));
            ac_within_ac_cz(ind2,nnn/10)=sum((distz<=nnn));
      
       end
end
%%

%%
% figure; plot((10:10:sdis),mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1),'g')
% hold on; plot((10:10:sdis),mean(ac_within_ac_c(:,:),1)./mean(meth_within_ac_c(:,:),1),'r')
% legend('methylation:acetylation contact ratio', 'acetylation:methylation contact ratio')
% xlabel('radius (nm)')
% ylabel('contact ratio')


% figure; plot((10:10:sdis),(mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1))./(max((mean(meth_withinc(:,:),1)./mean(ac_withinc(:,:),1)))),'g')

%%

%%
end