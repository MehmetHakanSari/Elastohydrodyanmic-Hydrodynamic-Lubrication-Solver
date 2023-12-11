function [K,Dlength,npuntiD,GFcut] = ElasticKer(E,nu1,L,npunti,dx,magn)
a = dx;
%definition of the elastic green function
fattore_incrementoGF = 2048; %2048
npuntiGF = npunti * fattore_incrementoGF;
dxGF = L / npuntiGF;
% narmtot=npuntiGF/2;
q0 = 2*pi / L;
four(1:npuntiGF) = 0;

for i = 0:(npuntiGF-1)
    
    if i > npuntiGF/2
        k = -npuntiGF + i;
    else
        k = i;
    end
    
    if k == 0
        four(i+1) = 0;
    else
        four(i+1)=((-1)^k)*(q0/(2*pi))*(2*sin((1/2)*k*q0*a)/(k*q0))*(-2*(1-nu1^2)*(1/abs(k*q0)));
        % four(i+1)=((-1)^k)*(1/abs(k*q0));
    end
end

Q = fft(four);
Q(npuntiGF+1) = Q(1);


x(1:npuntiGF+1) = 0; 

for i=1:npuntiGF+1
    x(i)=dxGF*(i-1-npuntiGF/2);
end


GFsovra(1,:) = x;
GFsovra(2,:) = Q/E;
% GFsovra(2,:)=Q;
GF(1:npunti) = 0;


for i = 1:npunti
    GF(1,i) = GFsovra(1, i + (fattore_incrementoGF-1) * (i-1));
    GF(2,i) = GFsovra(2, i + (fattore_incrementoGF-1) * (i-1));
end
GF(1,npunti+1) = GFsovra(1, npuntiGF+1);
GF(2,npunti+1) = GFsovra(2, npuntiGF+1);


Dlength = L / magn;     %1/4 di L
Rlength = Dlength / L;
npuntiD = round(npunti * Rlength);
% x2=linspace(-Dlength/2,Dlength/2,npuntiD);
GFcut = real(GF(2, round(length(GF) / 2) : round(length(GF)/2) + floor(npuntiD)-1));
% for i=1:npuntiD+1
% K(i,:) = real(GF(2,round(length(GF)/2)-i+1:round(length(GF)/2)+npuntiD-i+1));
% end
K = toeplitz(GFcut);
size(K)
end

