%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge Point Set Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract boundary points for the first shape
[shape_1_x,shape_1_y,shape_1_theta]=extractBoundary(shape_1_matrix);

if length(shape_1_x) >= nsamp
    [shape_1_x,shape_1_y,shape_1_theta] = sampleBoundaryPoints(shape_1_x,shape_1_y,shape_1_theta,nsamp); %Needs change
else  
    error('Shape 1: Insufficient samples')
end

%contour_1 is the 100 sample point matrix from shape 1 
contour_1 = [shape_1_x shape_1_y];

%Extract boundary points for the second shape
[shape_2_x,shape_2_y,shape_2_theta] = extractBoundary(shape_2_matrix);

if length(shape_2_x) >= nsamp
    [shape_2_x,shape_2_y,shape_2_theta] = sampleBoundaryPoints(shape_2_x,shape_2_y,shape_2_theta,nsamp);%Needs change
else
    error('Shape 2: Insufficient samples')
end

%contour_2 is the 100 sample point matrix from shape 2 
contour_2 = [shape_2_x shape_2_y]; 

%Display the contours

figure(1)
subplot(1,2,1)
imagesc(shape_1_matrix);
hold on
plot(contour_1(:,1),contour_1(:,2),'g^')
quiver(contour_1(:,1),contour_1(:,2),cos(shape_1_theta),sin(shape_1_theta),0.5,'g.')
hold off
axis('ij');
axis([1 shape_dim_2 1 shape_dim_1])

subplot(1,2,2)
imagesc(shape_2_matrix);
axis('image')
hold on
plot(contour_2(:,1),contour_2(:,2),'ro')
quiver(contour_2(:,1),contour_2(:,2),cos(shape_2_theta),sin(shape_2_theta),0.5,'r.')
hold off
axis('ij');
axis([1 shape_dim_2 1 shape_dim_1])

colormap(color_map)
drawnow

disp('extracting boundary points');

[x,y]=meshgrid(linspace(1,shape_dim_2,36),linspace(1,shape_dim_1,36));
x=x(:);y=y(:);M=length(x);