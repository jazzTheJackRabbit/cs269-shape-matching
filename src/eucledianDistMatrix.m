function d = eucledianDistMatrix(x, y)
[ndata, dimx] = size(x);
[ncentres, dimc] = size(y);
d = (ones(ncentres, 1) * sum((x.^2)', 1))' + ones(ndata, 1) * sum((y.^2)',1) - 2.*(x*(y'));
end