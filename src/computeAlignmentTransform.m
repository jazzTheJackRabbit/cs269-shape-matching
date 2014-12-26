%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Thin Plate Splines
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% warp each coordinate
fx_aff=cx(n_good+1:n_good+3)'*[ones(1,nsamp); contour_1'];
d2=max(eucledianDistMatrix(X3b,contour_1),0);
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
d2=max(eucledianDistMatrix(X3b,Xtan),0);
U=d2.*log(d2+eps);
fx_wrp=cx(1:n_good)'*U;
fx=fx_aff+fx_wrp;
fy_aff=cy(n_good+1:n_good+3)'*[ones(1,nsamp); Xtan'];
fy_wrp=cy(1:n_good)'*U;
fy=fy_aff+fy_wrp;

Ztan=[fx; fy]';
shape_1_theta=atan2(Ztan(:,2)-Z(:,2),Ztan(:,1)-Z(:,1));


figure(4)
plot(Z(:,1),Z(:,2),'g^',contour_2(:,1),contour_2(:,2),'ro');
axis('ij')
axis([1 shape_dim_2 1 shape_dim_1])

% show warped coordinate grid
fx_aff=cx(n_good+1:n_good+3)'*[ones(1,M); x'; y'];
d2=eucledianDistMatrix(X3b,[x y]);
fx_wrp=cx(1:n_good)'*(d2.*log(d2+eps));
fx=fx_aff+fx_wrp;
fy_aff=cy(n_good+1:n_good+3)'*[ones(1,M); x'; y'];
fy_wrp=cy(1:n_good)'*(d2.*log(d2+eps));
fy=fy_aff+fy_wrp;
hold on
plot(fx,fy,'.','markersize',12)
hold off
drawnow