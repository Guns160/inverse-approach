function y = li2(x)
%LI2 Real dilogarithm on [0,1], using the reflection around x=0.5.
values=double(x); flat=values(:); out=zeros(size(flat));
low=(flat>0)&(flat<=0.5); high=(flat>0.5)&(flat<1);
out(low)=charmonium.li2Series(flat(low));
z=1-flat(high);
out(high)=pi^2/6-log(flat(high)).*log1p(-flat(high))-charmonium.li2Series(z);
y=reshape(out,size(values));
end
