function [out]=cmatr(la,lb,n);
rows=la;
columns=lb;
out=zeros(la,lb);
for i=1:rows
    if n+1-(i-1) <= columns & n+1-(i-1)>0
       out(i,n+1-(i-1))=1;
    end
end
