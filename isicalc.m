function [out]=isicalc(hn,hn2,osr)
%
% ISICALC(hn,hn2,osr) calculates the intersymbol interference of the mached filter
% h (works also to matrix where a row presents a impulse response of the 
% transmit filter).
%  hn  =  transmit filter impulse response coeffs
%  hn2 =  receive filter impulse response coeffs
%  osr =  oversampleratio
% Initially written by Marko Kosunen 17.7.1998
%
% Last modified by: Marko Kosunen 28.8.1998

con=conv(hn,hn2);
l=size(con);
l=l(1,2);
b=(l-1)/2+1;

kmax=floor((b-1)/osr);
k1=-kmax:1:-1;
k2= 1:1:kmax;
k=[k1 k2].*osr;
 help3=abs(con(1,b));
   help2=sqrt(sum((abs(con(1,(k+b))./help3).^2)));
   help3=abs(con(1,b));
   out=20*log10(help2);
