import numpy as np
import scipy.io
import scipy.misc
import matplotlib.pyplot as plt
from PIL import Image

matInput = scipy.io.loadmat('digit_100_train_easy')
train_data = matInput['train_data']
label_train = matInput['label_train']

# % choose which two digits to compare:
mm=25;
nn=26;

# %%%Define flags and parameters:
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
annealing_rate=1; # annealing rate
w=4;
scale_factor=2.5; #scale factor

#Set the color map to gray 
plt.gray()

#Load the image
V1 = np.reshape(train_data[mm,:],(28,28)); #Reshape the row vector of image 1x784 into 28x28
V1 = scipy.misc.imresize(V1,scale_factor,'bilinear'); #Scale the image to 70x70 using bilinear interpolation

V2 = np.reshape(train_data[nn,:],(28,28));#Reshape the row vector of image 1x784 into 28x28
V2 = scipy.misc.imresize(V2,scale_factor,'bilinear'); #Scale the image to 70x70 using bilinear interpolation
[N1,N2] = V1.shape;#n1=70 , N2=70

if display_flag:
    #Display first image on a 2x2 grid, at position 1
    plt.figure(1)
    plt.subplot(2,2,1)
    plt.imshow(V1)
             
#   axis('image')
    plt.title(str(mm))

#   Display second image on a 2x2 grid, at position 2
    plt.figure(1)
    plt.subplot(2,2,2)
    plt.imshow(V2)
    plt.title(str(nn))
    
    plt.show()
#     axis('image')

