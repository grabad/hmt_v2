function [nucleus]=sim_nucleus_3d(nNucleoli)
px=10;
zrng=2;
zrng2=round(1000/px);

fact=75;
clear r;
N= 250;                                  % Number Of Sides To Polygon
a = sort(rand(N,1))*2*pi;
r(1)=3*rand+5;
for n=2:N
    r(n)=r(n-1)+0.04*randn;
    while r(n)<3
         r(n)=r(n-1)+0.04*randn;
    end
end
r=r';
xfact=(randn/2)+1;
yfact=(randn/2)+1;
while(yfact/xfact)>2||(yfact/xfact)<.5
    xfact=(randn/2)+1;
    yfact=(randn/2)+1;
end
while(yfact*xfact)<.25
    xfact=(randn/2)+1;
    yfact=(randn/2)+1;
end
for n102=1:zrng

x(:,n102) = xfact*cos(a*1).*r;
y(:,n102) = yfact*sin(a*1).*r;
x(:,n102)=x(:,n102)+abs(min(x(:,n102)));
y(:,n102)=y(:,n102)+abs(min(y(:,n102)));
while ((sum((round(x(:,n102)*fact)+50)>1400)))||((sum((round(y(:,n102)*fact)+50)>1400)))
    fact=fact-5;
end
rold=r;
dif=.1*(rand-.5);
r(1)=rold(1)+dif;
for n=2:N
    rold=r;
    dif=dif+.1*randn;
    r(n)=rold(n)+dif;
    while r(n)<3
         dif=dif+.1*randn;
         r(n)=rold(n)+dif;
         
    end
end
end

xtemp=zeros(N,zrng2);
ytemp=zeros(N,zrng2);
xtemp(:,1)=x(:,1);
ytemp(:,1)=y(:,1);
for n=1:length(x(:,1))
    ydif=y(n,1)-y(n,2);
    xdif=x(n,1)-x(n,2);
        for n2=2:zrng2

            xtemp(n,n2)=x(n,1)-(n2/zrng2)*xdif;
            ytemp(n,n2)=y(n,1)-(n2/zrng2)*ydif;
        end

end
x=xtemp;
y=ytemp;
%figure; plot([x; x(1)], [y; y(1)])
%%
nucleus=zeros(1400,1400,zrng2);


for n102=1:zrng2
    ntemp=zeros(1400,1400);
for n=1:N
    idxx(n)=round(x(n,n102)*fact)+50;
    idxy(n)=round(y(n,n102)*fact)+50;
    ntemp(idxx(n),idxy(n))=1;
end

while sum(size(ntemp)>1400)>0
    ntemp=zeros(1400);
    fact=fact-5;
    for n=1:N
        idxx(n)=round(x(n,n102)*fact)+50;
        idxy(n)=round(y(n,n102)*fact)+50;
        ntemp(idxx(n),idxy(n))=1;
    end

end
nucleus(:,:,n102)=ntemp;
s = size(nucleus(:,:,n102));
cn = round(s(1:2)/2); % [ycenter xcenter]
% find locations of outline pixels

for n=2:length(idxy)
    ydif=idxy(n)-idxy(n-1);
    xdif=idxx(n)-idxx(n-1);
    if max(abs(xdif),abs(ydif))>1
        for n2=1:(max(abs(xdif),abs(ydif))-1)
            co1=[idxx(n-1),idxy(n-1)];
            co2=[idxx(n),idxy(n)];
            xnew=idxx(n-1)+round(n2*(xdif/max(abs(xdif),abs(ydif))));
            ynew=idxy(n-1)+round(n2*(ydif/max(abs(xdif),abs(ydif))));
            nucleus(xnew,ynew,n102)=1;
        end

    end
end


    ydif=idxy(1)-idxy(n);
    xdif=idxx(1)-idxx(n);
    if max(abs(xdif),abs(ydif))>1
        for n2=1:(max(abs(xdif),abs(ydif))-1)
            co1=[idxx(1),idxy(1)];
            co2=[idxx(n),idxy(n)];
            xnew=idxx(n)+round(n2*(xdif/max(abs(xdif),abs(ydif))));
            ynew=idxy(n)+round(n2*(ydif/max(abs(xdif),abs(ydif))));
            nucleus(xnew,ynew,n102)=1;
          
        end
   
    end

nucleus(:,:,n102)=imfill(nucleus(:,:,n102));
se = strel('disk',5);
nucleus(:,:,n102)=imdilate(nucleus(:,:,n102),se);
end

%%
while 1
nucleolus=zeros(1400,1400,zrng2);
clear x; clear y;
for n99=1:nNucleoli
N= 75;                                  % Number Of Sides To Polygon
a = sort(rand(N,1))*2*pi;
clear r2
r2(1)=0.6*rand+1;
for n=2:N
    r2(n)=r2(n-1)+0.028*randn;
    while r2(n)<.3
         r2(n)=r2(n-1)+0.028*randn;
    end
end


r2=r2';
xfact=(randn/2)+1;
yfact=(randn/2)+1;

while(yfact/xfact)>2.5||(yfact/xfact)<.4
    xfact=(randn/2)+1;
    yfact=(randn/2)+1;
end
for n102=1:zrng
x(:,n102) = xfact*cos(a*1).*r2;
y(:,n102) = yfact*sin(a*1).*r2;

    
x(:,n102)=x(:,n102)+abs(min(x(:,n102)));
y(:,n102)=y(:,n102)+abs(min(y(:,n102)));
rold=r2;
dif=.04*(rand-.5);
r2(1)=rold(1)+dif;
for n=2:N
    rold=r2;
    dif=dif+.04*randn;
    r2(n)=rold(n)+dif;

end
end
xtemp=zeros(N,zrng2);
ytemp=zeros(N,zrng2);
xtemp(:,1)=x(:,1);
ytemp(:,1)=y(:,1);
for n=1:length(x(:,1))
    ydif=y(n,1)-y(n,2);
    xdif=x(n,1)-x(n,2);
        for n2=2:zrng2

            xtemp(n,n2)=x(n,1)-(n2/zrng2)*xdif;
            ytemp(n,n2)=y(n,1)-(n2/zrng2)*ydif;
        end

end
x=xtemp;
y=ytemp;
%figure; plot([x; x(1)], [y; y(1)])
%%

clear idxx
clear idxy

list=find(nucleus);
cent_ind=randi(length(list));
cent_nuc=list(cent_ind);
[I,J,K]=ind2sub([size(nucleus,1),size(nucleus,2),size(nucleus,3)],cent_nuc);
for n102=1:zrng2
    oli_temp=zeros(1400,1400);
for n=1:N
    idxx(n)=round(x(n,n102)*50)+1+I;
    idxy(n)=round(y(n,n102)*50)+1+J;
    oli_temp(idxx(n),idxy(n))=1;
end
while sum(size(oli_temp)>1400)>0
    oli_temp=zeros(1400);
    for n=1:N
        if idxx(n)>1000
             idxx(n)=round(x(n,n102)*50)-50;
        end
        if idxy(n)>1000
            idxy(n)=round(y(n,n102)*50)-50;
        end
         if idxx(n)<1
             idxx(n)=round(x(n,n102)*50)+50;
        end
        if idxy(n)<1
            idxy(n)=round(y(n,n102)*50)+50;
        end
        oli_temp(idxx(n),idxy(n))=1;
    end
end
nucleolus(:,:,n102)=nucleolus(:,:,n102)+oli_temp;
s = size(nucleolus);
cn = round(s(1:2)/2); % [ycenter xcenter]
% find locations of outline pixels

%[~,i]=min(pdist2([idxy,idxx],[idxy,idxx]),[],1)

for n=2:length(idxy)
    ydif=idxy(n)-idxy(n-1);
    xdif=idxx(n)-idxx(n-1);
    if max(abs(xdif),abs(ydif))>1
        for n2=1:(max(abs(xdif),abs(ydif))-1)
            co1=[idxx(n-1),idxy(n-1)];
            co2=[idxx(n),idxy(n)];
            xnew=idxx(n-1)+round(n2*(xdif/max(abs(xdif),abs(ydif))));
            ynew=idxy(n-1)+round(n2*(ydif/max(abs(xdif),abs(ydif))));
            nucleolus(xnew,ynew,n102)=1;
        end

    end
end


    ydif=idxy(1)-idxy(n);
    xdif=idxx(1)-idxx(n);
    if max(abs(xdif),abs(ydif))>1
        for n2=1:(max(abs(xdif),abs(ydif))-1)
            co1=[idxx(1),idxy(1)];
            co2=[idxx(n),idxy(n)];
            xnew=idxx(n)+round(n2*(xdif/max(abs(xdif),abs(ydif))));
            ynew=idxy(n)+round(n2*(ydif/max(abs(xdif),abs(ydif))));
            nucleolus(xnew,ynew,n102)=1;
          
        end
   
    end

nucleolus(:,:,n102)=imfill(nucleolus(:,:,n102));
end

end
if (sum(sum(nucleolus(:,:,1).*(1-nucleus(:,:,1))))<1) && (sum(sum(nucleolus(:,:,end).*(1-nucleus(:,:,end))))<1)
    break;
end

end

%%
nucleus=nucleus-nucleolus;
nucleus(nucleus<0)=0;
end