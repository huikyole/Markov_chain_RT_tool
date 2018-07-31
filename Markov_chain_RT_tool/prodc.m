function [c]=prodc(a,b)
%c=0;
%if a>=0
    c=prod([a+1:1:b+1])/(b+1);
%end
return