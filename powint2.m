function [out]=powint2(sig,pa,pl,na,nl)
% 
% POWINT(sig,n,lb) 
% Calculates the integrated noise power/signalpower ratio of the signal sig 
% 
%
% sig     = signal to be examined
% pa      = starpoint of the pass band
% pl      = stop poin  of the pass band
% na      = start of the stop band
% nl      = stop of the stop band
% 
%
%     
% Initially written by Marko Kosunen 17.7.1998
%
% Last modified by: Marko Kosunen 28.8.1998

a=abs(fft(sig));
d=sqrt(sum(a(1,pa:pl).^2)/(sum(a(1,na:nl).^2)));
e=-20*log10(d);
out=e;
