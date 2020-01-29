  function Cvals = pwCEval(a,b,c,d,x,z)
% 
% Pre:
%    a,b,c,d   column (n-1)-vectors that represent a piecewise
%              cubic polynomial.
%    x         column n-vector with x(1) < ... < x(n).
%    z         column m-vector with each z(i) in [x(1),x(n)].
%
% Post:
%    Cvals     column m-vector with the property that
%              Cvals(j) = C(z(j)), j=1:m. If z(j) is in 
%              the interval [x(i),x(i+1)], then 
%              C(z(j)) = a(i) + 
%                        b(i)(z(j)-x(i)) + 
%                        c(i)(z(j)-x(i))^2 +
%                        d(i)(z(j)-x(i))^3
% 
  m = length(z); 
  Cvals = zeros(m,1); 
  g=1;
  for j=1:m
     i = Locate(x,z(j),g);
     Cvals(j) = d(i)*(z(j)-x(i)) + c(i);
     Cvals(j) = Cvals(j)*(z(j)-x(i)) + b(i);
     Cvals(j) = Cvals(j)*(z(j)-x(i)) + a(i);
     g = i;
  end
 
