function out=csd2dec(csd)
s=size(csd);
l=s(1,2);
h=s(1,1);
out=zeros(1,h);
for i=1:h
  for j=1:l 
   out(1,i)=out(1,i)+csd(i,j)*1/2^(j-1);
  end
end

