function [blinking_data,gt]=blinking_sim_unrealistic(xpos,ypos,zpos,px,NA,nframes,size,identity2,fnameout)
%tub=tub(2:end-1,2:end-1);
load('count_dist');
load('fwhm');
for n=1:length(y647)
    cum647(n)=sum(y647(1:n));
    cum680(n)=sum(y680(1:n));
end
for n=1:length(Nf647)
    cumf647(n)=100*sum(Nf647(1:n));
    cumf680(n)=100*sum(Nf680(1:n));
end
sizea=(size*1000)/px; %k1=p3; k3=p3; k2=p2; k4=p6
frame_rate=10; %ms
%sizea=sizea(1);
%%
%notemp=nframes/nit;
%%
siz2=1000/px;
SRzoom=sizea/siz2;
sz=round(sizea/(px/SRzoom));
[x,y] = meshgrid(-sizea/2:sizea/2-1,-sizea/2:sizea/2-1);
xman=xpos;
yman=ypos;
spec_hetr=randn(length(xpos),1);
%PSF = genPSFparam(sz,px/SRzoom,NA,lambda)
gtx=[];
gty=[];
fkames=[];
speer=[];
int=[];
intesnity=[];
state=2*ones(length(xpos),1);
state=1;
frame_fin=[];
x_fin=[];
y_fin=[];
z_fin=[];
intens_fin=[];
centroid_fin=[]
for ind=1:length(xpos)
   randb=rand;
   if identity2(ind)>1.5 % 2 is 647, 1 is 680
    nblinks=find(cum647>=randb);
    nblinks=nblinks(1);
    for ind2=1:nblinks
       randf=rand;
       randi=rand;
       f_temp=find((randf)<=(cumf647));
       if isempty(f_temp)
           f_temp=0;
       end
       f_temp=(f_temp(1)+rand)*100;
       frame647(ind2)=round(f_temp);
       int_temp=find((randi)<=cumint647);
       int647(ind2)=int_temp(1);
       x_temp(ind2)=xpos(ind);
       y_temp(ind2)=ypos(ind);
       z_temp(ind2)=zpos(ind);
       centroid647(ind2)=671;
    end
        intens_fin=[intens_fin,int647];
        frame_fin=[frame_fin,frame647];
        x_fin=[x_fin,x_temp];
        y_fin=[y_fin,y_temp];
        z_fin=[z_fin,z_temp];
        centroid_fin=[centroid_fin,centroid647];
        clear x_temp; clear y_temp; clear z_temp; clear int647; clear frame647; clear centroid647;
   else
    nblinks=find(cum680>=randb);
    nblinks=nblinks(1);
        for ind2=1:nblinks
            randf=rand;
            randi=rand;
            f_temp=find((randf)<=(cumf680));
                if isempty(f_temp)
                    f_temp=0;
                end
            f_temp=(f_temp(1)+rand)*100;
            frame680(ind2)=round(f_temp);
            int_temp=find((randi)<=cumint680);
            int680(ind2)=int_temp(1);
            x_temp(ind2)=xpos(ind);
            y_temp(ind2)=ypos(ind);
            z_temp(ind2)=zpos(ind);
            centroid680(ind2)=705;
        end
        intens_fin=[intens_fin,int680];
        frame_fin=[frame_fin,frame680];
        x_fin=[x_fin,x_temp];
        y_fin=[y_fin,y_temp];
        z_fin=[z_fin,z_temp];
        centroid_fin=[centroid_fin,centroid680];
        
        clear x_temp; clear y_temp; clear z_temp; clear int680; clear frame680; clear centroid680;
   end
end
for ind3=1:length(z_fin)
   sig1(ind3)=zeroth_simulated(round(z_fin(ind3)+4100))/2.75;
   sig2(ind3)=first_simulated(round(z_fin(ind3)+4100))/2.75;

end
msig1=min(sig1);
msig2=min(sig2);
frame_fin=frame_fin-min(frame_fin)+1;
x_fin=x_fin-mean(x_fin)+1;
y_fin=y_fin-mean(y_fin)+1;
for ind=1:nframes
%    if mod((ind-1),notemp)==0
%       xpos=xman+rand-.5;
%       ypos=yman+rand-.5;
%    end

   on=find(frame_fin==ind);
   tempx=x_fin(on);
   tempy=y_fin(on);
   tempz=z_fin(on);
   tempf=frame_fin(on);
   tempi=intens_fin(on);
   tempc=centroid_fin(on);
   sig1t=sig1(on);
   sig2t=sig2(on);
   N=length(tempx);
   tub_blur=zeros(sizea);

   if mod(ind,100)==0
       perc=(ind/nframes)*100
   end

   for ind2=1:N
        xshift=tempx(ind2);
        yshift=tempy(ind2);
        %xt=x-xshift;%add x shift
        %yt=y-yshift;%add y shift
        Zo = sqrt(x.^2+y.^2);
        scale = sizea*px; 
        k_r = Zo./scale;
        Freq_max = NA/(tempc(ind2)*(sig1t(ind2)/msig1));           % cutting-off frequency
        pupil = k_r < Freq_max;         % pupil function
        PSFA = fftshift(fft2(pupil));               % fourier transform of the pupil
        PSF = PSFA.*conj(PSFA);
        h = PSF./sum(PSF(:));
        
        h=imtranslate(h,[xshift/px,yshift/px]);
        
        tub_blur=tub_blur+h*tempi(ind2);
   end

   
   
   blinking_data(:,:,ind)=(tub_blur);

   
   
end
gt=[(1:length(frame_fin))',frame_fin',x_fin',y_fin',z_fin',centroid_fin',sig1',sig2',intens_fin'];

%saveastiff(blinking_data, fnameout)
end