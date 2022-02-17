function sem = calcSEM(data,dim)
    if nargin == 1 && isvector(data)
        n = length(data);
        dim = 2;
    elseif nargin == 1 && ismatrix(data)
        n = size(data,1);
        dim = 1;
    else
        n = size(data,dim);
    end
    s = std(data,0,dim,'omitnan');
    sem = s/sqrt(n);
end