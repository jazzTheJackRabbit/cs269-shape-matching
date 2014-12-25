% calculate shape context cost
[a1,b1]=min(costmat,[],1);
[a2,b2]=min(costmat,[],2);
sc_cost=max(mean(a1),mean(a2));