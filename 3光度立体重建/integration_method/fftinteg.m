function z=fftinteg(p,q, cutoff)
% function z=fftinteg(p,q)
%
% computes a surface Z whose gradient (pz,qz):
%
%          ( dz  dz ) 
% (pz,qz)= ( --, -- )
%          ( dx  dy )
%
% is closest to p,q in a least-squares sense: 
%
% SUM_allpixels { (p-pz)^2 + (q-qz)^2 } --> minimum
% 
% Z,P,Q are M-by-N matrices. 
%
% Based on:
% Frankot, R.T and Chelappa, R.: "A Method for Enforcing Integrability 
% in Shape from Shading Algorithms", T-PAMI 10(4): 439-45, 1988. 
%
% (c) Ondrej Drbohlav, CMP Prague, 2003

P=fft2(p);
Q=fft2(q);
[m,n]=size(P);
[k,l]=meshgrid(0:n-1,0:m-1);
k(1,1)=1;
l(1,1)=1;
m2=floor(m/2);
n2=floor(n/2);

k(k>n2)=k(k>n2)-n;
l(l>m2)=l(l>m2)-m;

S=-i/(2*pi)*((k/n).*P+(l/m).*Q)./((k/n).^2+(l/m).^2);

%high-pass filtering of resulted bump-map JF
sS = fftshift(S);
pix = cutoff; 
r0 = round(size(S,1)/2);
c0 = round(size(S,2)/2);
for ii=r0-pix:r0+pix
  for jj=c0-pix:c0+pix
    if(sqrt((ii-r0)*(ii-r0)+(jj-c0)*(jj-c0)) < pix)
      sS(ii,jj) = 0;
    end;
  end;
end;
S = ifftshift(sS);

z=real(ifft2(S));
