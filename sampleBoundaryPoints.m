function [xi,yi,ti]=sampleBoundaryPoints(x,y,t,nsamp)
N=length(x);
k=3;
Nstart=min(k*nsamp,N);

randIdxs=randperm(N);
randIdxs=randIdxs(1:Nstart);

xi=x(randIdxs);
yi=y(randIdxs);
ti=t(randIdxs);


distMatrix=eucledianDistMatrix([xi yi],[xi yi]);
distMatrix=distMatrix+diag(Inf*ones(Nstart,1));

s=1;
while s
   % find closest pair
   [minVal1,rowIdxOfColMins]=min(distMatrix);
   [minVal2,minIdx]=min(minVal1);

   % remove one of the points
   xi(minIdx)=[];
   yi(minIdx)=[];
   ti(minIdx)=[];
   distMatrix(:,minIdx)=[];
   distMatrix(minIdx,:)=[];
   if size(distMatrix,1)==nsamp
      s=0;
   end
end

      