function idx = localMaxima(values)
%LOCALMAXIMA Interior weak-left/strict-right local maxima.
values=values(:); idx=find(values(2:end-1)>=values(1:end-2)& ...
    values(2:end-1)>values(3:end))+1;
end
