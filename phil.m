clear all; close all;
n = 1e2;
iseed = 34523;
rand('state',iseed);
spins = 1000;
hit = 0;
for m = 1:spins
u = rand(1,n);
usrt = sort(u,'ascend');
usrt = [0 usrt 1];
d = 3;
idx = 1:(n+2-d);
iu =  (1+d):(n+2);
S = usrt(iu) - usrt(idx);
M = max(S);
Mx(m)  = n*M - log(n) - (d-1)*log(log(n)) +log(factorial(d-1));
if (Mx(m) < 0)
hit = hit + 1;
end;
end;
ecdf(Mx);
xxx = -2:0.01:5;
G = exp(-exp(-xxx));
hold on;
plot(xxx,G,'r-.');
n = 1e5;
iseed = 27865;
rand('state',iseed);
spins = 1000;
hit = 0;
for m = 1:spins
u = rand(1,n);
usrt = sort(u,'ascend');
usrt = [0 usrt 1];
d = 3;
idx = 1:(n+2-d);
iu =  (1+d):(n+2);
S = usrt(iu) - usrt(idx);
M = max(S);
Mx(m)  = n*M - log(n) - (d-1)*log(log(n)) +log(factorial(d-1));
if (Mx(m) < 0)
hit = hit + 1;
end;
end;
hold on; 
ecdf(Mx);
grid on;
legend('n=100', 'Gumbel', 'n=100000');
