%Script to plot ShapeContexts
nPts = 4;
Pts = [1 34 72 80];
sp =1;
for idx =  1:nPts
b1 = BH1(Pts(idx),:);
b2 = BH2(cvec(Pts(idx)),:);
sc1 = zeros([12 5]);
sc2 = zeros([12 5]);
for k = 1:60
   sc1(sub2ind([12 5],k)) = b1(k); 
   sc2(sub2ind([12 5],k)) = b2(k); 
end

sc1n = uint8(sc1 .* (255/max(max(sc1))));
sc2n = uint8(sc2 .* (255/max(max(sc2)))); 
subplot(2,nPts,sp);imshow(sc1n); 
subplot(2,nPts,sp+nPts);imshow(sc2n);
sp = sp+1;
end