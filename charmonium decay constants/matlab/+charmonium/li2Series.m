function y = li2Series(z)
%LI2SERIES Power series for real 0<z<=0.5.
term=double(z(:)); total=term;
for k=2:10000
    term=term.*z(:); add=term/(k*k); updated=total+add;
    if all(abs(add) <= 4*eps(max(1,abs(updated))))
        y=updated; return;
    end
    total=updated;
end
error('Li2 power series did not converge');
end
