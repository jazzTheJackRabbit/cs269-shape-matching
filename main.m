
load parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load image from the Dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[train_data label_data] = loadDataset('./Dataset/train-images-idx3-ubyte','./Dataset/train-labels-idx1-ubyte',100,0);

% Choose label of the shape to compare:
digit_label_to_match  = randi([0, 9]);

digit_label_mask = label_data == digit_label_to_match;
index_of_label_to_match = find(digit_label_mask);

shape_1 = index_of_label_to_match(randi([1, length(index_of_label_to_match)]));

%Fix to not pick the same shape twice.
digit_label_mask(shape_1) = 0;
index_of_label_to_match = find(digit_label_mask);

shape_2 = index_of_label_to_match(randi([1, length(index_of_label_to_match)]));

color_map=flipud(gray);

[shape_1_matrix,shape_2_matrix,shape_dim_1,shape_dim_2] = loadShapes(train_data(:,:,shape_1),train_data(:,:,shape_2),sf);

edgeDetection

runShapeContextTPS

