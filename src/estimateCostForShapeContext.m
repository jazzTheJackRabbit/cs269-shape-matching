costmat_shape = computeHistogramCost(shapeContextHistogram1,shapeContextHistogram2);
theta_diff=repmat(shape_1_theta,1,nsamp)-repmat(shape_2_theta',nsamp,1);    
costmat_theta=0.5*(1-cos(theta_diff));    
costmat=(1-ori_weight)*costmat_shape+ori_weight*costmat_theta;
nptsd=nsamp+ndum;
costmat2=eps_dum*ones(nptsd,nptsd);
costmat2(1:nsamp,1:nsamp)=costmat;
cvec=hungarian(costmat2);