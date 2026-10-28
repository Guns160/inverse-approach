function test_pp_anomaly()
%TEST_PP_ANOMALY Check the Ward conversion and public PP option propagation.
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(fullfile(fileparts(here),'eta_c_pseudoscalar','matlab'));
addpath(fullfile(fileparts(here),'eta_c_axial_vector','matlab'));
mass=[2.9841;3.6378]; R=[1,-1,4;4,0,1]; mc=1.1;
[f,raw,shift,meta]=charmonium.ppDecayConstants(R,mass,mc);
assert(abs(shift(1)+0.079/mass(1)^2)<1e-15&&shift(2)==0);
assert(isequal(meta.applied,[true;false]));
assert(all(isnan(f(R<=0))));
valid=R>0; recovered=(f.*mass.^2-[-0.079;0])/(2*mc);
assert(max(abs(recovered(valid).^2-R(valid)))<1e-13);
[omitted,rawOmitted]=charmonium.ppDecayConstants(R,mass,mc,[],[]);
assert(isequaln(omitted,rawOmitted)&&isequaln(raw,rawOmitted));
[custom,~,customShift,customMeta]=charmonium.ppDecayConstants(R,mass,mc,0,0.02);
assert(all(customMeta.applied)&&customShift(1)==0);
assert(abs(customShift(2)-0.02/mass(2)^2)<1e-15);
assert(all(isnan(custom(R<=0))));
base=solve_eta_c_pseudoscalar(5,8,16,'DS',0.5);
other=solve_eta_c_pseudoscalar(5,8,16,'DS',0.5,'Anomaly1S',[],'Anomaly2S',0.02);
assert(isequaln(base.two_pole.residues_all,other.two_pole.residues_all));
assert(isequaln(base.two_pole.rho_all,other.two_pole.rho_all));
assert(base.two_pole.alpha_star==other.two_pole.alpha_star);
assert(isequal(other.metadata.anomaly.applied,[false;true]));
assert(isequaln(other.two_pole.f_all, ...
    other.two_pole.f_raw_all+[0;0.02/mass(2)^2]));
rejected=false;
try, solve_eta_c_axial_vector(5,8,16,'Anomaly1S',-0.079);
catch, rejected=true; end
assert(rejected,'AA must reject PP-only anomaly options');
disp('PP anomaly: Ward identity, raw values, custom inputs, unchanged inversion PASS');
end
