function r=nancorrcoef(x)

[~,N]=size(x);
r = zeros(N);
for i=1:N
    for j=i:N
        indx=find(~isnan(x(:,i)+x(:,j)));
        xi=x(indx,i)-mean(x(indx,i));
        xj=x(indx,j)-mean(x(indx,j));
        r(i,j) = (xi'*xj)./sqrt(xi'*xi)./sqrt(xj'*xj);
    end
end
r = r+r'-eye(N);
