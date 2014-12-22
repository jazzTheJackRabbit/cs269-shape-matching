% load in the digit database 
if ~(exist('train_data')&exist('label_train'))
load digit_100_train_easy;
%   load digit_100_train_hard;
end

% choose two digits to compare:
mm=35;
nn=42;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define flags and parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_flag=1;
affine_start_flag=1;
polarity_flag=1;
nsamp=100;
eps_dum=0.25;
ndum_frac=0.25;        
mean_dist_global=[];
ori_weight=0.1;
nbins_theta=12;
nbins_r=5;
r_inner=1/8;
r_outer=2;
tan_eps=1.0;
n_iter=6;
beta_init=1;
r=1; % annealing rate
w=4;
sf=2.5; %scale factor

cmap=flipud(gray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Image Loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[shape_1_matrix,shape_2_matrix,shape_dim_1,shape_dim_2] = load_shapes(mm,nn,train_data,sf)

if display_flag
    %Display first image on a 2x2 grid, at position 1
    figure(1)
    subplot(2,2,1)
    imagesc(shape_1_matrix);axis('image')
    title(int2str(mm))

    %Display second image on a 2x2 grid, at position 2
    subplot(2,2,3)
    imagesc(shape_2_matrix);axis('image')
    title(int2str(nn))

    %Show image in grayscale
    colormap(cmap)
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('extracting boundary points...')

%Extract boundary points for the first shape
[shape_1_x,shape_1_y,shape_1_theta]=boundary_extraction(shape_1_matrix);

nsamp1=length(shape_1_x);

if nsamp1 >= nsamp
    [shape_1_x,shape_1_y,shape_1_theta] = get_samples_1(shape_1_x,shape_1_y,shape_1_theta,nsamp); %Needs change
else  
    error('Shape 1: Insufficient samples')
end

%contour_1 is the 100 sample point matrix from shape 1 
contour_1 = [shape_1_x shape_1_y];

%Extract boundary points for the second shape
[shape_2_x,shape_2_y,shape_2_theta] = boundary_extraction(shape_2_matrix);

nsamp2 = length(shape_2_x);

if nsamp2 >= nsamp
    [shape_2_x,shape_2_y,shape_2_theta] = get_samples_1(shape_2_x,shape_2_y,shape_2_theta,nsamp);%Needs change
else
    error('Shape 2: Insufficient samples')
end

%contour_2 is the 100 sample point matrix from shape 2 
contour_2 = [shape_2_x shape_2_y];


%Display the contours
if display_flag
    subplot(2,2,2)
    plot(contour_1(:,1),contour_1(:,2),'g^')
    hold on
    quiver(contour_1(:,1),contour_1(:,2),cos(shape_1_theta),sin(shape_1_theta),0.5,'g.')
    hold off
    axis('ij');axis([1 shape_dim_2 1 shape_dim_1])
    title([int2str(length(shape_1_x)) ' samples'])
    subplot(2,2,4)
    plot(contour_2(:,1),contour_2(:,2),'ro')
    hold on
    quiver(contour_2(:,1),contour_2(:,2),cos(shape_2_theta),sin(shape_2_theta),0.5,'r.')
    hold off
    axis('ij');axis([1 shape_dim_2 1 shape_dim_1])
    title([int2str(length(shape_2_x)) ' samples'])
    drawnow	
end

if display_flag
    [x,y]=meshgrid(linspace(1,shape_dim_2,36),linspace(1,shape_dim_1,36));
    x=x(:);y=y(:);M=length(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute correspondences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%t1 is the angle of the gradient for image 1
k=1;
s=1;
ndum=round(ndum_frac*nsamp);
out_vec_1=zeros(1,nsamp);
out_vec_2=zeros(1,nsamp);

%Compute correspondence and alignment transform for each iteration
while s
    disp(['iter=' int2str(k)])
    disp('computing shape contexts for (deformed) model...')
    [BH1,mean_dist_1]=sc_compute(contour_1',zeros(1,nsamp),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec_1);
    disp('done.')
    
    % apply the scale estimate from the warped model to the test shape
    disp('computing shape contexts for target...')
    [BH2,mean_dist_2]=sc_compute(contour_2',zeros(1,nsamp),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec_2);
    disp('done.')

    if affine_start_flag
        if k==1
            % use huge regularization to get affine behavior
            lambda_o=1000;
        else
            lambda_o=beta_init*r^(k-2);	 
        end
    else
        lambda_o=beta_init*r^(k-1);
    end
    beta_k=(mean_dist_2^2)*lambda_o;

    costmat_shape=hist_cost_2(BH1,BH2);
    theta_diff=repmat(shape_1_theta,1,nsamp)-repmat(shape_2_theta',nsamp,1);
    %   costmat_theta=abs(atan2(sin(theta_diff),cos(theta_diff)))/pi;
    if polarity_flag
        % use edge polarity
        costmat_theta=0.5*(1-cos(theta_diff));
    else
        % ignore edge polarity
        costmat_theta=0.5*(1-cos(2*theta_diff));
    end      
    costmat=(1-ori_weight)*costmat_shape+ori_weight*costmat_theta;
    nptsd=nsamp+ndum;
    costmat2=eps_dum*ones(nptsd,nptsd);
    costmat2(1:nsamp,1:nsamp)=costmat;
    cvec=hungarian(costmat2);
    %   cvec=hungarian_fast(costmat2);

    % update outlier indicator vectors
    [a,cvec2]=sort(cvec);
    out_vec_1=cvec2(1:nsamp)>nsamp;
    out_vec_2=cvec(1:nsamp)>nsamp;

    X2=NaN*ones(nptsd,2);
    X2(1:nsamp,:)=contour_1;
    X2=X2(cvec,:);
    X2b=NaN*ones(nptsd,2);
    X2b(1:nsamp,:)=contour_1;
    X2b=X2b(cvec,:);
    Y2=NaN*ones(nptsd,2);
    Y2(1:nsamp,:)=contour_2;

    % extract coordinates of non-dummy correspondences and use them
    % to estimate transformation
    ind_good=find(~isnan(X2b(1:nsamp,1)));
    n_good=length(ind_good);
    X3b=X2b(ind_good,:);
    Y3=Y2(ind_good,:);

    if display_flag
        figure(2)
        plot(X2(:,1),X2(:,2),'g^',Y2(:,1),Y2(:,2),'ro')
        hold on
        h=plot([X2(:,1) Y2(:,1)]',[X2(:,2) Y2(:,2)]','k-');

        if 1
            %	 set(h,'linewidth',1)
            quiver(contour_1(:,1),contour_1(:,2),cos(shape_1_theta),sin(shape_1_theta),0.5,'g')
            quiver(contour_2(:,1),contour_2(:,2),cos(shape_2_theta),sin(shape_2_theta),0.5,'r')
        end
        
        hold off
        axis('ij')
        title([int2str(n_good) ' correspondences (warped X)'])
        axis([1 shape_dim_2 1 shape_dim_1])
        drawnow	
    end

    if display_flag
        % show the correspondences between the untransformed images
        figure(3)
        plot(contour_1(:,1),contour_1(:,2),'g^',contour_2(:,1),contour_2(:,2),'ro')
        ind=cvec(ind_good);
        hold on
        plot([X2b(:,1) Y2(:,1)]',[X2b(:,2) Y2(:,2)]','k-')
        hold off
        axis('ij')
        title([int2str(n_good) ' correspondences (unwarped X)'])
        axis([1 shape_dim_2 1 shape_dim_1])
        drawnow	
    end

    [cx,cy,E]=bookstein(X3b,Y3,beta_k);

    % calculate affine cost
    A=[cx(n_good+2:n_good+3,:) cy(n_good+2:n_good+3,:)];
    s=svd(A);
    aff_cost=log(s(1)/s(2));

    % calculate shape context cost
    [a1,b1]=min(costmat,[],1);
    [a2,b2]=min(costmat,[],2);
    sc_cost=max(mean(a1),mean(a2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Thin Plate Splines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % warp each coordinate
    fx_aff=cx(n_good+1:n_good+3)'*[ones(1,nsamp); contour_1'];
    d2=max(dist2(X3b,contour_1),0);
    U=d2.*log(d2+eps);
    fx_wrp=cx(1:n_good)'*U;
    fx=fx_aff+fx_wrp;
    fy_aff=cy(n_good+1:n_good+3)'*[ones(1,nsamp); contour_1'];
    fy_wrp=cy(1:n_good)'*U;
    fy=fy_aff+fy_wrp;

    Z=[fx; fy]';

    % apply the warp to the tangent vectors to get the new angles
    Xtan=contour_1+tan_eps*[cos(shape_1_theta) sin(shape_1_theta)];
    fx_aff=cx(n_good+1:n_good+3)'*[ones(1,nsamp); Xtan'];
    d2=max(dist2(X3b,Xtan),0);
    U=d2.*log(d2+eps);
    fx_wrp=cx(1:n_good)'*U;
    fx=fx_aff+fx_wrp;
    fy_aff=cy(n_good+1:n_good+3)'*[ones(1,nsamp); Xtan'];
    fy_wrp=cy(1:n_good)'*U;
    fy=fy_aff+fy_wrp;

    Ztan=[fx; fy]';
    shape_1_theta=atan2(Ztan(:,2)-Z(:,2),Ztan(:,1)-Z(:,1));

    if display_flag
        figure(4)
        plot(Z(:,1),Z(:,2),'g^',contour_2(:,1),contour_2(:,2),'ro');
        axis('ij')
        title(['k=' int2str(k) ', \lambda_o=' num2str(lambda_o) ', I_f=' num2str(E) ', aff.cost=' num2str(aff_cost) ', SC cost=' num2str(sc_cost)])
        axis([1 shape_dim_2 1 shape_dim_1])
        
        % show warped coordinate grid
        fx_aff=cx(n_good+1:n_good+3)'*[ones(1,M); x'; y'];
        d2=dist2(X3b,[x y]);
        fx_wrp=cx(1:n_good)'*(d2.*log(d2+eps));
        fx=fx_aff+fx_wrp;
        fy_aff=cy(n_good+1:n_good+3)'*[ones(1,M); x'; y'];
        fy_wrp=cy(1:n_good)'*(d2.*log(d2+eps));
        fy=fy_aff+fy_wrp;
        hold on
        plot(fx,fy,'k.','markersize',1)
        hold off
        drawnow
    end

    % update contour_1 for the next iteration
    contour_1=Z;

    if k==n_iter
        s=0;
    else
        k=k+1;
    end
end

%%%
%%% grayscale warping
%%%

% [x,y]=meshgrid(1:N2,1:N1);
% x=x(:);y=y(:);M=length(x);
% fx_aff=cx(n_good+1:n_good+3)'*[ones(1,M); x'; y'];
% d2=dist2(X3b,[x y]);
% fx_wrp=cx(1:n_good)'*(d2.*log(d2+eps));
% fx=fx_aff+fx_wrp;
% fy_aff=cy(n_good+1:n_good+3)'*[ones(1,M); x'; y'];
% fy_wrp=cy(1:n_good)'*(d2.*log(d2+eps));
% fy=fy_aff+fy_wrp;
% disp('computing warped image...')
% V1w=griddata(reshape(fx,N1,N2),reshape(fy,N1,N2),V1,reshape(x,N1,N2),reshape(y,N1,N2));   
% fz=find(isnan(V1w)); V1w(fz)=0;
% ssd=(V2-V1w).^2;
% ssd_global=sum(ssd(:));
% if display_flag
%    figure(5)
%    subplot(2,2,1)
%    im(V1)
%    subplot(2,2,2)
%    im(V2)
%    subplot(2,2,4)
%    im(V1w)
%    title('V1 after warping')
%    subplot(2,2,3)
%    im(V2-V1w)
%    h=title(['SSD=' num2str(ssd_global)]);
%    colormap(cmap)
% end
% 
% %%%
% %%% windowed SSD comparison
% %%%
% wd=2*w+1;
% win_fun=gaussker(wd);
% % extract sets of blocks around each coordinate
% % first do 1st shape; need to use transformed coords.
% win_list_1=zeros(nsamp,wd^2);
% for qq=1:nsamp
%    row_qq=round(contour_1(qq,2));
%    col_qq=round(contour_1(qq,1));
%    row_qq=max(w+1,min(N1-w,row_qq));
%    col_qq=max(w+1,min(N2-w,col_qq));
%    tmp=V1w(row_qq-w:row_qq+w,col_qq-w:col_qq+w);
%    tmp=win_fun.*tmp;
%    win_list_1(qq,:)=tmp(:)';
% end
% % now do 2nd shape
% win_list_2=zeros(nsamp,wd^2);
% for qq=1:nsamp
%    row_qq=round(Y(qq,2));
%    col_qq=round(Y(qq,1));
%    row_qq=max(w+1,min(N1-w,row_qq));
%    col_qq=max(w+1,min(N2-w,col_qq));
%    tmp=V2(row_qq-w:row_qq+w,col_qq-w:col_qq+w);
%    tmp=win_fun.*tmp;
%    win_list_2(qq,:)=tmp(:)';
% end
% ssd_all=sqrt(dist2(win_list_1,win_list_2));
% if 0
%    % visualize paired-up patches
%    for qq=1:nsamp
%       im([reshape(win_list_1(qq,:),wd,wd) reshape(win_list_2(b2(qq),:),wd,wd); ...
% 	  reshape(win_list_2(qq,:),wd,wd) reshape(win_list_1(b1(qq),:),wd,wd)])
%       colormap(cmap)
%       pause
%    end
% end
% % loop over nearest neighbors in both directions, project in
% % both directions, take maximum
% cost_1=0;
% cost_2=0;
% for qq=1:nsamp
%    cost_1=cost_1+ssd_all(qq,b2(qq));
%    cost_2=cost_2+ssd_all(b1(qq),qq);
% end
% ssd_local=(1/nsamp)*max(mean(cost_1),mean(cost_2));
% ssd_local_avg=(1/nsamp)*0.5*(mean(cost_1)+mean(cost_2));
% if display_flag
% %   set(h,'string',['local SSD=' num2str(ssd_local) ', avg. local SSD=' num2str(ssd_local_avg)])
%    set(h,'string',['local SSD=' num2str(ssd_local)])
% end

