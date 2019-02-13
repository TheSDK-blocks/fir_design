a=rrcos(0.22,2,37);
b=a.*(kaiser(37,4))';
load ntclagrideals
c=chfil;

ideal=rrcos(0.22,2,1001);
%c=lagrfdes(ideal,37,0.5,0.61,2,10^(75/10),10^(45/10));
g=splot([zpad([a;b;c],1024);zpad(ideal,1024)],[-90 0]);

as=g(1,1:512);
bs=g(2,1:512);
cs=g(3,1:512);

ideals=g(4,1:512);

isia=isicalc(a,ideal,2)
isib=isicalc(b,ideal,2)
isic=isicalc(c,ideal,2)
%isic2=isicalc(c2,ideal,2)
isiideal=isicalc(ideal,ideal,2)

aclra=powint2(zpad(a,1024),1,256,313,513)
aclrb=powint2(zpad(b,1024),1,256,313,513)
aclrc=powint2(zpad(c,1024),1,256,313,513)
%aclrc2=powint2(zpad(c2,1024),1,256,313,513)
aclrideal=powint2(zpad(ideal,1024),1,256,313,513)

%plot section
clf
xmin=1;
xmax=512;
ymin=-120;
ymax=5;
handle=axes;
h2=gcf;
x=1:512;
f=plot(x,as,'k--',x,bs,'k:',x,cs,'k-',x,ideals,'-.k');
t4=legend(f,'Truncation','Window','Lagrange','1001 coeffs');
fs=14;
set(t4,'FontSize' , fs);
set(t4,'Linewidth',2)
axis([1 xmax ymin ymax])
grid on
set(handle,'XTick',[128 256 313  384])
set(handle,'XTicklabel',['0.25'; '0.5 '; '0.61'; '0.75'])
set(handle,'YTick',[-120:10:0])
set(handle,'Linewidth',2)
set(f,'Linewidth',2)
t1=title('Comparison of filter design methods');
set(t1,'FontSize' , fs);
t2=ylabel('Relative power [dB]');
set(t2,'FontSize' , fs);
t3=xlabel('Normalized frequency');
set(t3,'FontSize' , fs);
%print -deps filter_comparison.eps
%print -deps filter_comparison_2.eps
print -deps filter_comparison_3.eps
