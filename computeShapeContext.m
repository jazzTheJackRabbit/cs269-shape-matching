%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Shape Contexts
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeContextHistogram,mean_dist]=computeShapeContext(xy_samples,theta_samples,mean_dist,nbins_theta,nbins_r,r_bin_hist_start,r_bin_hist_end,out_vec);

nsamp=size(xy_samples,2);
in_vec = out_vec==0;

% Compute r,theta arrays
radial_distance_array = real(sqrt(eucledianDistMatrix(xy_samples',xy_samples')));                                          
theta_array_abs = atan2(xy_samples(2,:)'*ones(1,nsamp)-ones(nsamp,1)*xy_samples(2,:),xy_samples(1,:)'*ones(1,nsamp)-ones(nsamp,1)*xy_samples(1,:))';
theta_array = theta_array_abs-theta_samples'*ones(1,nsamp);

% Normalize distance by mean and ignore outliers
if isempty(mean_dist)
   tmp=radial_distance_array(in_vec,:);
   tmp=tmp(:,in_vec);
   mean_dist=mean(tmp(:));
end
r_array_normalized = radial_distance_array/mean_dist;

% Bin the distances on a log scale
r_bin_hist_edges = logspace(log10(r_bin_hist_start),log10(r_bin_hist_end),nbins_r);
hist_r_bin_index_array = zeros(nsamp,nsamp);

for m=1:nbins_r
   hist_r_bin_index_array = hist_r_bin_index_array+(r_array_normalized<r_bin_hist_edges(m));
end

flag_zero_array=hist_r_bin_index_array>0; % flag all points inside outer boundary

theta_array_2 = rem(rem(theta_array,2*pi)+2*pi,2*pi);
hist_theta_bin_index_array = 1+floor(theta_array_2/(2*pi/nbins_theta));

nbins=nbins_theta*nbins_r;
shapeContextHistogram=zeros(nsamp,nbins);

% Combine the theta and radial distance arrays into a sparse vector
for i=1:nsamp
   fzn=flag_zero_array(i,:)&in_vec;
   Sn=sparse(hist_theta_bin_index_array(i,fzn),hist_r_bin_index_array(i,fzn),1,nbins_theta,nbins_r);
   shapeContextHistogram(i,:)=Sn(:)';
end




