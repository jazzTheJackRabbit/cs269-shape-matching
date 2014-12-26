% calculate affine cost
A=[cx(n_good+2:n_good+3,:) cy(n_good+2:n_good+3,:)];
warping=svd(A);
aff_cost=log(warping(1)/warping(2));