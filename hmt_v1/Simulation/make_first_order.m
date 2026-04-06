%split=4;
%eta=1;
%Pbg = [1200];
%Ibgp = 1*sum(sum(mean(data,3)));
%EM = 100;
[n11,n11,Nf]=size(data);
CF680=csvread('CF 680.csv');
Alexa_Fluor_647=csvread('Alexa_Fluor_647.csv');
% spec647=spec647-min(spec647);
% spec647=spec647/sum(spec647)
% spec680=spec680-min(spec680);
% spec680=spec680/sum(spec680)
% zero_647=zero_647-min(zero_647);
% zero_647=zero_647/sum(zero_647)
% zero_680=zero_680-min(zero_680);
% zero_680=zero_680/sum(zero_680);
% blank=zeros(length(spec647));
% blank(ceil(length(spec647)/2),:)=spec647;
% blank2=zeros(length(zero_647));
% blank2(:,ceil(length(zero_647)/2))=zero_647;

spt647 = Alexa_Fluor_647(:,3);
w1 = Alexa_Fluor_647(:,1);
spt680 = CF680(:,3);
w2 = CF680(:,1);
n2 = length(spt647);
%figure; plot(w1,spt/max(spt),'r');%,wl,spt2/max(spt2));
load('new_cal.mat')
dif=900:1100;
centroid=p1.*dif.^2+p2.*dif+p3;
wl2 = centroid; % for Alexa 647 (99.8% signal)
n3 = length(wl2);
sptimg647 = interp1(w1,spt647,wl2)'; %Make signal correct for given spectral dispersion
sptimg647(isnan(sptimg647))=[];
sptimg647 = sptimg647/sum(sptimg647(:)); %makes signal intensity add up to 1
sptimg680 = interp1(w2,spt680,wl2)'; %Make signal correct for given spectral dispersion
sptimg680(isnan(sptimg680))=[];
sptimg680 = sptimg680/sum(sptimg680(:)); %makes signal intensity add up to 1

% first_ord=zeros(150);
% sz1=length(first_ord);
% [sz2x,sz2y]=size(sptimg2(:,:,1));
% first_ord((sz1/2-(sz2x/2-1)):(sz1/2+(sz2x/2)),(sz1/2-(sz2y/2-1)):(sz1/2+(sz2y/2)))=sptimg2(:,:,1);