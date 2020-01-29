% (b)
load shuo.mat

x = cell2mat(x(:));
y = cell2mat(y(:));

figure(1);
subplot(2,1,1);
plot(x,y,'-');
title('Piecewise Linear Plot');
axis(v);
xlabel('x-coordinate');
ylabel('y-coordinate');


N = length(x);
t =  zeros(N,1);

for i=2:length(x)
  t(i)= t(i-1) + sqrt((x(i)-x(i-1)).^2 + (y(i)-y(i-1)).^2);
end

[xa,xb,xc,xd] = mySpline(t, x);
[ya,yb,yc,yd] = mySpline(t, y);

fac=100;
tref = zeros(1, fac*(N-1) + 1);
for k = 1:N-1
 i = fac*(k-1)+1;
 dt = t(k+1)-t(k);
 for l = 0:fac-1
  tref(i+l) = t(k) + l*dt/fac;
 end
end
tref(fac*(N-1)+1) = t(N);

xx = pwCEval(xa, xb, xc, xd, t, tref); 
yy = pwCEval(ya, yb, yc, yd, t, tref);


subplot(2,1,2);
plot(xx,yy);
title('Cubic Spline Interpolation');
axis(v);
xlabel('x-coordinate');
ylabel('y-coordinate');


% (a)
function [a,b,c,d] = mySpline(t,x)
    n = length(x);
    T = zeros(n);
    s = zeros(n,1);
    r = zeros(n,1);
    
    a = zeros(n-1,1);
    b = zeros(n-1,1);
    c = zeros(n-1,1);
    d = zeros(n-1,1);
    
    dx = zeros(n-1,1);
    for i = 1:n-1
        dx(i) = t(i+1) - t(i);
    end
    
    y_ = zeros(n-1,1);
    for j = 1:n-1
        y_(j) = (x(j+1) - x(j))/dx(j);
    end
    
    T(1,1) = dx(2)^2;
    T(1,2) = dx(2)^2-dx(1)^2;
    T(1,3) = -dx(1)^2;
    T(n,n-2) = (dx(n-1))^2;
    T(n,n-1) = dx(n-1)^2 + dx(n-2)^2;
    T(n,n) = -dx(n-2)^2;
    
    r(1) = 2*(dx(2)).^2*y_(1)-2*(dx(1)).^2*y_(2);
    r(n) = 2*(dx(n-1)).^2*y_(n-2)-2*(dx(n-1)).^2*y_(n-2);
    
    for k=2:n-1
        T(k,k-1) = dx(k);
        T(k,k) = 2*(dx(k)+dx(k-1));
        T(k,k+1) = dx(k-1);
    end
    
    for l=2:n-1
        r(l) = 3*(dx(l)*y_(l-1)+dx(l-1)*y_(l));
    end
    
    s = T\r;
    
    for m=1:n-1
        a(m) = x(m);
        b(m) = s(m);
        c(m) = (3*y_(m)-2*s(m)-s(m+1))/dx(m);
        d(m) = (s(m+1)+s(m)-2*y_(m))/((dx(m))^2);
    end
end
