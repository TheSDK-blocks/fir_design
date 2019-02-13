%lagrfdes(transfil,ntaps,fpass,fstop,osr,bs,ci)
%Function to optimize the matched filter pair
% Maximizes the costfuction
% PP-bs*PS-ci*ISI where PP=Passband energy PS=stopband energy
% and ISI is intersymbol interference.
%
% transfil= the fixed filter of the macthed filter pair
% ntaps   = number of taps in filter to be designed (odd)
% fpass   = pass band corner frequency (relative to nygvist) 
% fstop   = stop band corner frequency (relative to nyqvist)
% osr     = over sample ratio=fsamp
% bs      = weight factor of the stop band energy
% ci      = weight factor of the inter symbol interference
%
% Initially written by: Marko Kosunen 
% Last modified by Marko Kosunen 19.1.1999

function [out]=lagrfdes(transfil,ntaps,fpass,fstop,osr,bs,ci);

%Some help parameters
s1=size(transfil);  
lt=s1(1,2);
transfil=transfil'; %Basic form of the filter is column vector
for i=1:ntaps+lt-1

   A(:,i)=cmatr(ntaps,lt,i-1)*transfil; %Convolution matrix.
%Multiplication with this matrix equals to convolution.

end
% calculation of the ISI matrix B. (Hr'*B)*(Hr'*B)'=ISI^2 
l=ntaps+lt-1;
b=(l-1)/2+1;
kmax=floor((b-1)/osr);
k1=-kmax:1:-1;
k2= 1:1:kmax;
k=[k1 k2].*osr;

B=A(:,k+b);
clear k;

AC=A(:,b) % Hr'*AC= Center Tap 

% Calculation of the stopband (PS) and passband (PP) energies.  
for i=1:ntaps
    for k=1:ntaps
       r=i-1;
       c=k-1;
       if r-c==0
         PP(i,k)=2*fpass;
         PS(i,k)=osr-2*fstop;
       else
          PP(i,k)=(osr/(pi*(r-c)))*sin(2*pi*fpass*(r-c)/osr);
          PS(i,k)=(osr/(pi*(r-c)))*(sin(pi*(r-c))-sin(2*pi*fstop*(r-c)/osr));
       end
    end
end

W=PP;
D=PS;
% Help function Q
Q=2*W-2*bs*D-2*ci*B*B';
QINV=inv(Q);

Hr=QINV*AC/((QINV*AC)'*AC); %Optimized filter.
lambda=-1/((QINV*AC)'*AC);
Fisi=(Hr'*B)';
F=(Hr'*A)';
EP=Hr'*W*Hr %Passband energy
ES=Hr'*D*Hr %Stopband energy
SNR=10*log10(EP/ES) %signal to noise ratio
test=Hr'*AC-1 %Just a test
ISI=10*log10(Fisi'*Fisi) %ISI
out=Hr'; % row vector Hr' is the designed filter
PS
