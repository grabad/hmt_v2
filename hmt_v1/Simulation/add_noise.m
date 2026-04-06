split=2;
eta=1;
%Pbg = [2000];

EM = 1;
[n11,n11,Nf]=size(blinking_first);
bgimgp = repmat(ones(n11,n11),[1 1 Nf])*Ibgp;
bgimg = repmat(ones(n11,n11),[1 1 Nf])*Ibg;
Nsp = random('Poisson',data);
Ns = random('Poisson',blinking_first);
Nbgp = random('Poisson',bgimgp);
Nbg = random('Poisson',bgimg);
Nr = random('norm', 0,RN,n11,n11,Nf);
Nrp = random('norm', 0,RN,n11,n11,Nf);
first_noise = sqrt(2-1/EM)*((Ns - blinking_first) + (Nbg - bgimg))+ Nr/EM + blinking_first;
data_noise = sqrt(2-1/EM)*((Nsp - data) + (Nbgp - bgimgp))+ Nrp/EM + data;

% imwrite(uint16(data_noise(:,:,1)),'noise_dat.tif')
% for ind=2:Nf
% imwrite(uint16(data_noise(:,:,ind)),'noise_dat.tif','WriteMode','append');
% end