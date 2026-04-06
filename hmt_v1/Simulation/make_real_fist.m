function [blinking_first,gtfist]=make_real_fist(first_ord,spec_het,gt,FOV,px,NA,nframes)
load('new_cal')
sizea=(FOV*1000)/px;
first_ord=nonzeros(first_ord);
frames=gt(:,2);
[frames,idx]=sort(frames);
gt=gt(idx,:);
frame=1;
blinking_first=zeros((FOV*1000)/px,(FOV*1000)/px,nframes);
[x,y] = meshgrid(-sizea/2:sizea/2-1,-sizea/2:sizea/2-1);
for ind=1:length(gt(:,1))
        while gt(ind,2)>frame
           frame=frame+1; 
           if frame==394
               ff=4
           end
        end
        if frame<=nframes
        int=gt(ind,9);
        peak=gt(ind,6);
        sig2t=gt(ind,8);
        dif=(sqrt(4*p1*(peak-p3)+p2^2)-p2)/(2*p1);
        xshift=gt(ind,3);
        yshift=gt(ind,4);
        %xt=x-xshift;%add x shift
        %yt=y-yshift;%add y shift
        msig1=min(gt(ind,7));
        Zo = sqrt(x.^2+y.^2);
        scale = sizea*px; 
        k_r = Zo./scale;
        Freq_max = NA/(gt(ind,6)*(sig2t/msig1));           % cutting-off frequency
        pupil = k_r < Freq_max;         % pupil function
        PSFA = fftshift(fft2(pupil));               % fourier transform of the pupil
        PSF = PSFA.*conj(PSFA);
        h = PSF./sum(PSF(:));
        
        h=imtranslate(h,[xshift/px+(dif-950),yshift/px]);
        sptimg = conv2(h,first_ord','same');
        blinking_first(:,:,frame)=blinking_first(:,:,frame)+sptimg*int;
      
        else
        end
end
% sizea=length(first_ord);
% gtx=[];
% gty=[];
% intesnity=[];
% for ind=1:max(unique(gt(:,3)))
%     if mod(ind,100)==0
%         (ind/(max(unique(gt(:,3)))))*100
%     end
%     temp1=gt(:,1);
%     temp2=gt(:,2);
%     temp3=gt(:,4);
%     temp4=gt(:,5);
%     gttemp=[temp1(gt(:,3)==ind),temp2(gt(:,3)==ind),temp3(gt(:,3)==ind),temp4(gt(:,3)==ind)];
%     rfirst=zeros(sizea);
%     clear shifth;
%     clear int;
%    for ind2=1:length(gttemp(:,1))
%        int(ind2)=3*gttemp(ind2,4);
%       
%        shifth(ind2)=spec_het/res*gttemp(ind2,3);
%        shiftx=gttemp(ind2,1)-sizea/2+shifth(ind2);
%        shifty=gttemp(ind2,2)-sizea/2;
%        h1=imtranslate(first_ord,[shiftx,shifty]);
%        rfirst=rfirst+h1*int(ind2);
%    end
%    
%    blinking_first(:,:,ind)=(rfirst);
%    if ~isempty(gttemp(:,1))
%        gtx=[gtx;gttemp(:,1)+shifth'];
%        gty=[gty;gttemp(:,2)];
%        intesnity=[intesnity;int'];
%        
%    
%    else
%        hal=1;
%    end
%    %fkames=[fkames;ind*ones(length(tempx),1)];
% end
% gtfist=[gtx,gty,intesnity];
% end