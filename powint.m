function [out]=powint(sig,n,lb)
% 
% POWINT(sig,n,lb) 
% Calculates the integrated noise power/signalpower ratio of the signal sig 
% 
%
% sig     = signal to be examined
% n       = terms to be taken to  account
% lb      = lower bound of signal band
% 
%
%     
% Initially written by Marko Kosunen 17.7.1998
%
% Last modified by: Marko Kosunen 28.8.1998

a=abs(fft(sig));
b=a(1,1:n);
c=sort(b);
d=sqrt(sum(c(1,1:n-lb).^2)/(sum(c(1,n-lb+1:n).^2)));
e=20*log10(d);
d=c./max(c);
out=e;