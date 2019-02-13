function [out]=csdc(a,nbits,mterms,mode)
% 
% CSDC(a,nbits,mterms,mode) 
% optimates the csd filter coefficients
% using Samueli optimization algorithm
%
% a       = original filter coefficients
% nbits   = number of nonzero bits in csd presentation
% 
% mterms  = number of terms in csd code
% mode    =  mode=1 canonic signed digit presentation, mode/=1 minimum signed
%            digit presentation (not canonic)
%     
% Initially written by Marko Kosunen 17.7.1998
%
% Last modified by: Marko Kosunen 28.8.1998



s=size(a);
l=s(1,2);
for i=1:l
  a_a=abs(a(1,i));
  siga=sign(a(1,i));
  h=zeros(1,mterms);
  ts=1;
  n=0;  
  
  if (mode == 2 | mode==1) & a_a > 0.5
  nbits2=nbits+1;
  else
  nbits2=nbits;
  end  

 
  if a_a > csd2dec(ones([1 nbits2]))
    fprintf(1, '%s %s', 'Domain' ,'overflow');
  return
  end
  
  while sum(abs(h)) < nbits2 & a_a > 1/2^mterms,
  
    while 1/2^n > a_a,
      n=n+1;     
    end
  
  if  a_a > 1
    h(n+1)=ts;
    a_a=abs(a_a-1/2^n);
  

  elseif abs(a_a-1/2^(n-1)) <= abs(a_a-1/2^(n))
    h(n)=ts;
    a_a=abs(a_a-1/2^(n-1));
    ts=-ts;
  
  elseif n < mterms
    h(n+1)=ts;
    a_a=abs(a_a-1/2^(n));
  
  
  
  elseif abs(a_a) < 1/2^mterms
  a_a=0;
  
  else
  
  h(n)=ts;
  a_a=0;
  end
 
 end
h=siga.*h;


out(i,:)=h;

end

if mode == 1
changes=1;

while changes == 1
changes=0;

for k=1:l

  if abs(csd2dec(out(k,:)))<=1


    for i=1:mterms-1  
      
      if abs(out(k,i))==1 & abs(out(k,i+1))==1
      changes=1;  
       out(k,i-1)=out(k,i);
       out(k,i)=0;
       out(k,i+1)=-out(k,i+1);
      end
    end
   end
end
end
end









