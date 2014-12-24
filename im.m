function im(X,vec);
if nargin==1
   imagesc(X)
else
   imagesc(X,vec)
end
pixval on
title(inputname(1))
colormap(gray)
colorbar
axis('image')

