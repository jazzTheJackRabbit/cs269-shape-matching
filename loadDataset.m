%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load Dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [images labels] = loadDataset(imageFile, labelFile, numberOfDigitsToRead, offset)
    
    % Read digits
    fid = fopen(imageFile, 'r', 'b');
    header = fread(fid, 1, 'int32');
    
    %Hardcoded Header value.
    if header ~= 2051
        error('Invalid image file header');
    end
    
    %Number of images in the dataset
    count = fread(fid, 1, 'int32');
    if count < numberOfDigitsToRead+offset
        error('Trying to read too many digits');
    end
    
    %Number of rows for each image: 32 bit integer
    numberOfRows = fread(fid, 1, 'int32');
    
    %Number of columns for each image: 32 bit integer
    numberOfColumns = fread(fid, 1, 'int32');
    
    %Offset to start from
    if offset > 0
        fseek(fid, numberOfColumns*numberOfRows*offset, 'cof');
    end
    
    %Init the return image array
    images = zeros([numberOfRows numberOfColumns numberOfDigitsToRead]);
    
    %Traverse the dataset file and read each image into return image array
    for i=1:numberOfDigitsToRead
        for y=1:numberOfRows
            images(y,:,i) = fread(fid, numberOfColumns, 'uint8');
        end
    end
    
    fclose(fid);

    % Read digit labels
    fid = fopen(labelFile, 'r', 'b');
    header = fread(fid, 1, 'int32');
    
    %Hardcoded Header value.
    if header ~= 2049
        error('Invalid label file header');
    end
    
    %Number of labels in the dataset
    count = fread(fid, 1, 'int32');
    if count < numberOfDigitsToRead+offset
        error('Trying to read too many digits');
    end
    
    %Offset to start from
    if offset > 0
        fseek(fid, offset, 'cof');
    end
    
    %Read readDigit number of labels from the label's file
    labels = fread(fid, numberOfDigitsToRead, 'uint8');
    fclose(fid);
    
     %Normalize the points so that they are anti aliased
     images = normalizePixelValue(images);
    
end

function digits = normalizePixelValue(digits)
    digits = double(digits);
    for i=1:size(digits, 3)
        digits(:,:,i) = digits(:,:,i)./255.0;
    end
end
