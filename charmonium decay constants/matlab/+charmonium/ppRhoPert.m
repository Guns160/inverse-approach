function rho = ppRhoPert(s,mc,alphaS)
%PPRHOPERT Sum of P-P LO and NLO spectral densities at supplied RG inputs.
[lo,nlo]=charmonium.ppRhoComponents(s,mc,alphaS); rho=lo+nlo;
end
