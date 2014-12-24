%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Resize shapes into 70x70
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shape_1,shape_2,shapes_size_dim1,shapes_size_dim2] = loadShapes(shape_1, shape_2,scaleFactor)
    shape_1=imresize(shape_1,scaleFactor,'bil'); %Scale the image to 70x70 using bilinear interpolation
    shape_2=imresize(shape_2,scaleFactor,'bil'); %Scale the image to 70x70 using bilinear interpolation
    [shapes_size_dim1,shapes_size_dim2]=size(shape_1);%n1=70 , N2=70
end