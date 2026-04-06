function [meth_pos,ac_pos]=calc_decay(xm,ym,zm,xa,ya,za)
[xm,order]=sort(xm); ym=ym(order); zm=zm(order);
[xa,order]=sort(xa); ya=ya(order); za=za(order);


meth_within=zeros(length(xm),25);
ac_within=zeros(length(xm),25);
sdis=1000;
for ind2=1:length(xm)
       mask=(abs(xm-xm(ind2))<=sdis).*(abs(ym-ym(ind2))<=sdis).*(abs(zm-zm(ind2))<=sdis);
       mask=logical(mask); 
           
       [ind,dist]=rangesearch([xm(mask),ym(mask),zm(mask)],[xm(ind2),ym(ind2),zm(ind2)],sdis);
        for nnn=10:10:sdis
            meth_within(ind2,nnn/10)=sum((dist{1}<=nnn))-1;
        end
        
        mask=(abs(xa-xm(ind2))<=sdis).*(abs(ya-ym(ind2))<=sdis).*(abs(za-zm(ind2))<=sdis);
        mask=logical(mask);    
       [ind,dist]=rangesearch([xa(mask),ya(mask),za(mask)],[xm(ind2),ym(ind2),zm(ind2)],sdis);
        for nnn=10:10:sdis
            ac_within(ind2,nnn/10)=sum((dist{1}<=nnn));
        end
end

meth_within_ac=zeros(length(xa),25);
ac_within_ac=zeros(length(xa),25);
for ind2=1:length(xa)
       mask=(abs(xm-xa(ind2))<=sdis).*(abs(ym-ya(ind2))<=sdis).*(abs(zm-za(ind2))<=sdis);
       mask=logical(mask);
  
       [ind,dist]=rangesearch([xm(mask),ym(mask),zm(mask)],[xa(ind2),ya(ind2),za(ind2)],sdis);
        for nnn=10:10:sdis
            meth_within_ac(ind2,nnn/10)=sum((dist{1}<=nnn));
        end
        
        mask=(abs(xa-xa(ind2))<=sdis).*(abs(ya-ya(ind2))<=sdis).*(abs(za-za(ind2))<=sdis);
        mask=logical(mask);
   
       [ind,dist]=rangesearch([xa(mask),ya(mask),za(mask)],[xa(ind2),ya(ind2),za(ind2)],sdis);
        for nnn=10:10:sdis
            ac_within_ac(ind2,nnn/10)=sum((dist{1}<=nnn))-1;
        end
end


sdis=1000;
pix_meth_in_meth=zeros(1,length(xm));
pix_meth_in_ac=zeros(1,length(xm));


%%
meth_ac_con=mean(meth_within(:,:),1)./mean(ac_within(:,:),1);
ac_meth_con=mean(ac_within_ac(:,:),1)./mean(meth_within_ac(:,:),1);
[meth_pos,ac_pos]=compile_function_in(ac_meth_con,meth_ac_con);


end