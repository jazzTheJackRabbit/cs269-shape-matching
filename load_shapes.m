function [V1,V2,N1,N2]=load_shapes(shape_1, shape_2,train_data,sf);
    V1=reshape(train_data(shape_1,:),28,28)'; %Reshape the row vector of image 1x784 into 28x28
    V1=imresize(V1,sf,'bil'); %Scale the image to 70x70 using bilinear interpolation
    V2=reshape(train_data(shape_2,:),28,28)';%Reshape the row vector of image 1x784 into 28x28
    V2=imresize(V2,sf,'bil'); %Scale the image to 70x70 using bilinear interpolation
    [N1,N2]=size(V1);%n1=70 , N2=70
end