function [x,y,theta,samplePoints]=extractBoundary(im);

samplePoints=contourc(im,[.5 .5]);
[Gx,Gy]=gradient(im);

samplePoints(:,samplePoints(1,:)==0.5)=[];
nSamples=size(samplePoints,2);

theta=zeros(nSamples,1);

for n=1:nSamples
   x_rnd=round(samplePoints(1,n));
   y_rnd=round(samplePoints(2,n));
   theta(n)=atan2(Gy(y_rnd,x_rnd),Gx(y_rnd,x_rnd))+pi/2;
end

x=samplePoints(1,:)';
y=samplePoints(2,:)';

end