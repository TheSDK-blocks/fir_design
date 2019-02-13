function [out]=csdopt(hn,nbits,csdwidth,omegap,omegas,omegas2,alpha,osr,isir,isibound,powbound,nstep,scalstep,mode)
% 
% CSDOPT(hn,nbits,csdwidth,omegap,omegas,omegas2,alpha,osr,isir,isibound,powbound,scalstep,mode) 
% optimates the csd filter coefficients
% using Samueli optimization algorithm
%
% hn      = original filter coefficients
% nbits   = number of nonzero bits in csd presentation
% 
% alpha   = roll-off factor
% omegap  = pass band cut-off frequency
% omegas  = stop band start frequency
% osr     = oversample ratio
% isir    = intersymbol interference to stopband attenuation ratio 
%             (att=isir*isi)
% isibound= boundary for ISI in winbdow methods
% powbound= boundary for stop band attenuation in window methods
% scalstep= minimum step used in scaling
% mode    = optimization mode mode=1 uses samueli optimization mode=2
%           minimizes the time domain squared error sum.
%           mode3 and mode 5 are window methods 
%           CAUTION: MODE 2 DOES NOT USE BIVARIATE SEARCH
% Initially written by Marko Kosunen 17.7.1998
%
% Last modified by: Marko Kosunen 30.12.1998

hnopt=hn;
s=size(hn);
l=s(1,2);

%scalstep=2^(-csdwidth-1);

ideal=rrcos(alpha,osr,1001);

if mode==3

 a1=0;
 a3=5;
 a2=(a1+a3)/2;
 win1=kaiser(l,a1)';
 win2=kaiser(l,a2)';
 win3=kaiser(l,a3)';
 hnoptv=hn;

 for i=1:nstep 
   win1=kaiser(l,a1)';
   win2=kaiser(l,a2)';
   win3=kaiser(l,a3)';
   hnopt1=hn.*win1;
   hnopt2=hn.*win2;
   hnopt3=hn.*win3;
   hnopt11=csd2dec(csdc(hnopt1,nbits,csdwidth,2));
   hnopt22=csd2dec(csdc(hnopt2,nbits,csdwidth,2));
   hnopt33=csd2dec(csdc(hnopt3,nbits,csdwidth,2));
   isi1=isicalc(hnopt11,ideal,osr)
   isi2=isicalc(hnopt22,ideal,osr)
   isi3=isicalc(hnopt33,ideal,osr) 
   sppowb=powint2(zpad(hnoptv,1024),1,omegap,omegas,omegas2)

   if isi2 <= isibound & isi3 >= isibound
      a1=a2
      a2=(a2+a3)/2
      hnopt=hnopt2;
      hnoptv=hnopt22;
   elseif isi1 <=isibound & isi2 >= isibound
      a3=a2
      a2=(a1+a3)/2
      hnopt=hnopt1;
      hnoptv=hnopt11;
   elseif isi3 < isibound
      a3=a3+1
      a2=(a1+a3)/2;
      hnopt=hnopt3;
      hnoptv=hnopt33;
      sppowb=powint2(zpad(hnoptv,1024),1,omegap,omegas,omegas2)
   end
 end
end

if mode == 5

 a1=0;
 a3=5;
 a2=(a1+a3)/2;
 hnopt=hn;
 hnoptv=hn;
 sppowb=powint2(zpad(hnopt,1024),1,omegap,omegas,omegas2)

 for i=1:nstep 
    win1=kaiser(l,a1)';
    win2=kaiser(l,a2)';
    win3=kaiser(l,a3)';
    hnopt1=hn.*win1;
    hnopt2=hn.*win2;
    hnopt3=hn.*win3;
    hnopt11=csd2dec(csdc(hnopt1,nbits,csdwidth,2));
    hnopt22=csd2dec(csdc(hnopt2,nbits,csdwidth,2));
    hnopt33=csd2dec(csdc(hnopt3,nbits,csdwidth,2));
    pow1=powint2(zpad(hnopt11,1024),1,omegap,omegas,omegas2)
    pow2=powint2(zpad(hnopt22,1024),1,omegap,omegas,omegas2)
    pow3=powint2(zpad(hnopt33,1024),1,omegap,omegas,omegas2)

    if pow2 >= powbound & pow3 <= powbound
      a1=a2
      a2=(a2+a3)/2
      hnopt=hnopt3;
      hnoptv=hnopt33;
    elseif pow1 >=powbound & pow2 <= powbound
      a3=a2
      a2=(a1+a3)/2
      hnopt=hnopt1;
      hnoptv=hnopt11;
    elseif pow3 < powbound
        return
    end
 end
end
if mode==3 | mode ==5
a=isicalc(hnoptv,ideal,osr)
b=powint2(zpad(hnoptv,1024),1,omegap,omegas,omegas2)
end
if a > isibound | b > powbound
    return
end

if mode==1

  sigmaprev=0;
  hnopt=hn;

  for i=1:-scalstep:1/2+scalstep

    hscal=hnopt.*i;
    hscalcsd=csdc(hscal,nbits,csdwidth,2);
    hscalcsddec=csd2dec(hscalcsd);
    hscalzpad=zpad(hscalcsddec,1024);
    isi=isicalc(hscalcsddec,ideal,osr);
    hscalfreq=abs(fft(hscalzpad));
    gain=mean(hscalfreq(1,1:512*omegap));
    attn=max(20*log10(hscalfreq(1,omegas*512:512)/gain));
    sigma=max(attn,isi/isir);
  
    if sigma < sigmaprev

      hscalcsdopt=hscalcsd;
      sigmaprew=sigma;
    end

  end

hnopt=hscalcsdopt;

end

if mode == 2
 hscalcsdopt=csdc(hn,nbits,csdwidth,2); 
 helper=csd2dec(csdc(hn,nbits,csdwidth,2));
 errorprew=sum(abs(hn-helper).^2);

 for i=1:-scalstep:1/2+scalstep

  hscal=hn.*i;
  hscalcsd=csdc(hscal,nbits,csdwidth,2);
  hscalcsddec=csd2dec(hscalcsd);
  error=sum((abs(hscal-hscalcsddec)./i).^2);
  
  if error < errorprew
   errorprew=error;
   hscalcsdopt=hscalcsd;
 end

end

out=hscalcsdopt;
return
end



if mode==3 | mode==5

hscalcsdopt=csdc(hnopt,nbits,csdwidth,2);
hscalcsddec=csd2dec(hscalcsdopt);
sppow=b;

for i=1:-scalstep:1/2+scalstep
  hscal=hnopt.*i;
  hscalcsd=csdc(hscal,nbits,csdwidth,2);
  hscalcsddec=csd2dec(hscalcsd);
  sppowb=powint2(zpad(hscalcsddec,1024),1,omegap,omegas,omegas2);
  isib=isicalc(ideal,hscalcsddec,osr);
 
 if sppowb < sppow & isib <=isibound
   sppow=sppowb;
   hscalcsdopt=hscalcsd;
  end

end

hnopt=hscalcsdopt;

end

if mode ~=2
%bivariate search

hnopt=csd2dec(hnopt);
sppowb=powint2(zpad(hnopt,1024),1,omegap,omegas,omegas2);
hnoptb=hnopt;
change=1;
s=size(hnoptb);
l=s(1,2);
while change==1;

change=0;

for i=1:(l-1)/2
 for j=1:(l-1)/2+1;
  for k=1:2
   for m=1:2 
   
   if k==1       
   hnoptb(1,i)=hnoptb(1,i)-1/2^csdwidth;
   hnoptb=[hnoptb(1,1:(l-1)/2+1) fliplr(hnoptb(1,1:(l-1)/2))];
   elseif k==2
   hnoptb(1,i)=hnoptb(1,i)+1/2^csdwidth;
   hnoptb=[hnoptb(1,1:(l-1)/2+1) fliplr(hnoptb(1,1:(l-1)/2)) ];
   end
  
   if m==1 & i~=j       
   hnoptb(1,j)=hnoptb(1,j)-1/2^csdwidth;
   hnoptb=[hnoptb(1,1:(l-1)/2+1) fliplr(hnoptb(1,1:(l-1)/2))];
   elseif m==2
   hnoptb(1,j)=hnoptb(1,j)+1/2^csdwidth;
   hnoptb=[ hnoptb(1,1:(l-1)/2+1) fliplr(hnoptb(1,1:(l-1)/2)) ];
   end 
   
   hnoptb=csd2dec(csdc(hnoptb,nbits,csdwidth,2));
   sppow=powint2(zpad(hnoptb,1024),1,omegap,omegas,omegas2);
   isib=isicalc(hnoptb,ideal,osr);
   
   if sppow < sppowb & isib <= isibound
   sppowb = sppow;
   hnopt=hnoptb;
   change=1;
   end  
  
  end
 end
end
end
end
end
sppowb

isi=isicalc(hnopt,ideal,osr)
out=csdc(hnopt,nbits,csdwidth,2);
 
 
 






 


























































