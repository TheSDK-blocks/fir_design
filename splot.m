% USAGE splot(sig,yborders) plots the spectrum of the signal sig
% sig = matrix where row corresponds one signal to be plotted
% yborders= [YMIN YMAX]

function [out]=splot(sig,yborders)
a=size(sig);
b=a(1,1);
c=a(1,2);
d=floor(c/2);           
hmat(1:6,:)=str2mat('y-', 'y--', 'm-', 'm--', 'c-', 'c--'); 
hmat(7:14,:)=str2mat('b-', 'b--', 'g-', 'g--', 'w-', 'w--', 'r-', 'r--');
figure(1);
clf;
for i=1:b
 hold on;
 h(i,:)=fft(sig(i,:));
 h(i,:)=abs(h(i,:));
 h(i,:)=20*log10(h(i,:)/max(h(i,:)));
 l=find(isinf(h(i,:)));
 h(i,l)=yborders(1,1)*ones(size(l));
 plot([1:d],h(i,1:d),sprintf('%s',hmat((rem(i-1,14)+1),1)));
end
axis([1 d yborders]);
hold off;
out=h;
