function [gtx,gty,gtz,f1x,f1y,f1z,identity2]=sim_anti(tempx,tempy,zpos,identity)
ind2=1;
sdis=4.9; %Make sure no 2 H3 histones are within 6 nm to avoid overlapping
% for n=1:length(tempx)
% %     
%         mask=(abs(tempx-tempx(ind2))<=sdis).*(abs(tempy-tempy(ind2))<=sdis).*(abs(zpos-zpos(ind2))<=sdis);
%         mask=logical(mask);
%   if sum(mask)>1
%        [ind,dist]=knnsearch([tempx(mask),tempy(mask),zpos(mask)],[tempx(ind2),tempy(ind2),zpos(ind2)],'K',2);
%         if dist(2)<6
%            
%            tempx(ind2)=[]; tempy(ind2)=[];zpos(ind2)=[];
%            ind2=ind2-1;
%            
%         end
%   end
%         ind2=ind2+1;
%         
% end
gtx=tempx; gty=tempy; gtz=zpos;
efficiency=0.25;
anti_length=9;
lucky=zeros(1,length(gtx));
ind3=0;
ant_dist=sqrt(3^2+5^2);
for n=1:length(gtx)
    if rand<efficiency
        ind3=ind3+1;
        lucky(n)=1;
        identity2(ind3)=identity(n);
        theta1=2*pi()*rand-pi();
        phi1=(3*pi()/2)*rand-3*pi()/4;
        theta2=2*pi()*rand-pi();
        phi2=2*pi()*rand-pi();
        xin=rand;
        yin=rand;
        zin=rand;
        magn=sqrt(xin^2+yin^2+zin^2);
        xshift(ind3)=ant_dist*(xin/magn);
        yshift(ind3)=ant_dist*(yin/magn);
        zshift(ind3)=ant_dist*(zin/magn);
        ax(ind3)=gtx(n)+xshift(ind3)+anti_length*cos(theta1)*cos(phi1);
        ay(ind3)=gty(n)+yshift(ind3)+anti_length*sin(theta1)*cos(phi1);
        az(ind3)=gtz(n)+zshift(ind3)+anti_length*sin(phi1);
        f1x(ind3)=ax(ind3)+anti_length*cos(theta2)*cos(phi2);
        f1y(ind3)=ay(ind3)+anti_length*sin(theta2)*cos(phi2);
        f1z(ind3)=az(ind3)+anti_length*sin(phi2);
    end
end
inx=gtx(find(lucky));
iny=gty(find(lucky));
inz=gtz(find(lucky));
deletion=zeros(1,length(f1x));
for n=1:length(f1x)
   anticorx=(inx(n)+xshift(n)):((ax(n)-inx(n)-xshift(n))/10):ax(n);
   anticory=(iny(n)+yshift(n)):((ay(n)-iny(n)-yshift(n))/10):ay(n);
   anticorz=(inz(n)+zshift(n)):((az(n)-inz(n)-zshift(n))/10):az(n);
   anticorx=[anticorx(1:10),ax(n):((f1x(n)-ax(n))/10):f1x(n)];
   anticory=[anticory(1:10),ay(n):((f1y(n)-ay(n))/10):f1y(n)];
   anticorz=[anticorz(1:10),az(n):((f1z(n)-az(n))/10):f1z(n)];
   gtxt=gtx(1:end ~= n); gtyt=gty(1:end ~= n); gtzt=gtz(1:end ~= n);
    mask=(sqrt(abs(gtxt-anticorx).^2+abs(gtyt-anticory).^2+abs(gtzt-anticorz).^2)<=sdis);
    mask=logical(mask);
    if sum(sum(mask))>0
        deletion(n)=1;
    end
end
f1x(find(deletion))=[];
f1y(find(deletion))=[];
f1z(find(deletion))=[];


end

