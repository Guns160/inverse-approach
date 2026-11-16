function [shape,gamma,partial] = pp_bw_spectrum(s,model)
%PP_BW_SPECTRUM Running BW in s, without finite-area renormalization.
% shape has units GeV^-2; gamma and partial widths have units GeV.
phi=pp_bw_phase_space(s,model);
partial=phi./model.phi_pole.*model.pole_partial_widths;
gamma=model.rest_width+sum(partial,2);
z=model.mass*gamma;
shape=z./(pi*((double(s(:))-model.mass^2).^2+z.^2));
end
