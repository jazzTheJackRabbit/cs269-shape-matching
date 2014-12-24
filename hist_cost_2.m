function histogramCost = computeHistogramCost(BH1,BH2);
% HC=hist_cost_2(BH1,BH2);
%
% same as hist_cost.m but BH1 and BH2 can be of different lengths

[nsamp1,nbins]=size(BH1);
[nsamp2,nbins]=size(BH2);

BH1n=BH1./repmat(sum(BH1,2)+eps,[1 nbins]); %BH1 ./ sumBH1, where sumBH1 is the repeated matrix of BH1 summed along the rows and for each element added with the eps and this vector is repeated 60 times.
BH2n=BH2./repmat(sum(BH2,2)+eps,[1 nbins]);
tmp1=repmat(permute(BH1n,[1 3 2]),[1 nsamp2 1]);
tmp2=repmat(permute(BH2n',[3 2 1]),[nsamp1 1 1]);
histogramCost=0.5*sum(((tmp1-tmp2).^2)./(tmp1+tmp2+eps),3);
