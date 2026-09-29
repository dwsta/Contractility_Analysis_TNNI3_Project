function [EArr,nuArr,GArr,KArr,LinE,NH] = func_runElastographyVideos(X,Y,UB,VB,inds,runLinE,runNH)

% Wrapper for elastography of cardiomyocyte videos

% RS: Assume deformations are already smooth
%
% %% Filter high frequency with FFT and also find the derivatives
% % divergence of displacment dv
% dVB = nan(size(UB)); 
% % Filtered displacements
% UBF = nan(size(UB)); VBF = nan(size(UB)); 
% 
% for i = 1 : size(UB,3)
% 
%     Utt = inpaint_nans(UB(:,:,i),3);
%     Vtt = inpaint_nans(VB(:,:,i),3);
%     [Ufilt,Vfilt,dUdx,dVdy] = FourierFiltAndDiff(X,Y,Utt,Vtt);
%     UBF(:,:,i) = Ufilt; VBF(:,:,i) = Vfilt; 
%     dVB(:,:,i) = dUdx + dVdy;
%     
% end
UBF = UB; VBF = VB;
dVB = divergence_rik(X,Y,UB,VB);
%% Fit (E,nu) and Fourier tractions


% inds = [40 90 120]; % Indices of the peak in the time series
% Postinds = [6 22 38 55 71 87 104 120];

if nargin <= 5
    runNH = 1;
    runLinE = 1;
end

[LinE, NH] = runoptimNHLinE(UBF,VBF,UBF,VBF,dVB,X,Y,inds,runLinE,runNH);

EArr=[];nuArr=[];
GArr=[];KArr=[];

if runLinE
for i = 1: length(inds)
EArr(i)=LinE(i).sol.E;
nuArr(i)=LinE(i).sol.nu;
end; end

if runNH
for i = 1: length(inds)
GArr(i)=NH(i).sol.G;
KArr(i)=NH(i).sol.K;
end; end