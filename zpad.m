% USAGE zpad(sig,npoints) zero pads the signal in matrix sig. (row extension)
% sig     = signal matrix
% npoints = total signal points after padding

function [out]=zpad(sig,npoints)
a=size(sig);
b=a(1,1);
c=a(1,2);
out=zeros([b npoints]);
for i=1:b
  out(i,1:c)=sig(i,:);
end
 
