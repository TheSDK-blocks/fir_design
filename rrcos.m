function [out]=rrcos(alpha,m,n)
%
% RRCOS - Function to calculate the inpulse response of
%         root raised cosine filter
% rrcos(alpha,m,n)
% alpha = roll-off factor 
% m     = oversample ratio
% n     = number of filter coefficients (odd)
% initially written by Marko Kosunen
% last modifified 27.6.1998
% by Marko Kosunen



b=1;

for i=-(n-1)/2:(n-1)/2

  t=i/m;

  if t==0

   h=1;

   elseif abs(t)==1/(4*alpha)

    term1=4*alpha/(pi*(1-alpha)+4*alpha);

    term4=(pi-2)*cos(pi/(4*alpha));

    term5=(pi+2)*sin(pi/(4*alpha));

    h=term1*(term4+term5)/(4*sqrt(2));

  else

   term1=4*alpha/(pi*(1-alpha)+4*alpha);

   term2=(4*alpha*t)*(1-(4*alpha*t)^2);

   term3=sin(pi*(1-alpha)*t);

   term4=4*alpha*t*cos(pi*(1+alpha)*t);

   h=term1/term2*(term3+term4);

  end

out(1,b)=h;

b=b+1;

end


