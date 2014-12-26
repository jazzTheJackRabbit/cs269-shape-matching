%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Histogram Cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function histogramCost = computeHistogramCost(shapeContext1,shapeContext2);

[nsamp1,nbins] = size(shapeContext1);
[nsamp2,nbins] = size(shapeContext2);

shapeContext1_normalized = shapeContext1./repmat(sum(shapeContext1,2)+eps,[1 nbins]); 
shapeContext2_normalized = shapeContext2./repmat(sum(shapeContext2,2)+eps,[1 nbins]);
tmp1 = repmat(permute(shapeContext1_normalized,[1 3 2]),[1 nsamp2 1]);
tmp2 = repmat(permute(shapeContext2_normalized',[3 2 1]),[nsamp1 1 1]);

histogramCost=0.5*sum(((tmp1-tmp2).^2)./(tmp1+tmp2+eps),3);
