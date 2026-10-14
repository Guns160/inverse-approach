function mass = massForChannel(phys,name)
%MASSFORCHANNEL Look up the two physical pole masses for the channel.
for row=1:size(phys.mass_table,1)
    if strcmp(phys.mass_table{row,1},name), mass=phys.mass_table{row,2}; return; end
end
error('Unknown channel %s',name);
end
