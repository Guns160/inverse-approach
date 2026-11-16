function info = pp_bw_state_info(channel)
%PP_BW_STATE_INFO Input registry only; a reserved state is not a width model.
if nargin<1, channel='5S'; end
assert((ischar(channel)&&isrow(channel))||(isstring(channel)&&isscalar(channel)), ...
    'Channel must be 3S, 4S, 5S, or 6S');
channel=upper(char(channel));
assert(any(strcmp(channel,{'3S','4S','5S','6S'})), ...
    'Channel must be 3S, 4S, 5S, or 6S');
info=struct('channel',channel,'status','available','usage','','mass',NaN, ...
    'reference_mass',NaN,'truncation_note','','blocker','');
switch channel
    case '3S'
        info.mass=4.0119; info.reference_mass=4.007;
        info.usage='diagnostic: too close to the two delta poles; poor extraction';
        info.truncation_note='Only DDstar is modeled; the high-side DstarDstar channel is omitted';
    case '4S'
        info.mass=4.2693; info.reference_mass=4.276;
        info.usage='alternative';
        info.truncation_note='Four-channel truncation; DD0star coupling is unavailable';
    case '5S'
        info.mass=4.4508; info.reference_mass=4.47;
        info.usage='primary';
        info.truncation_note='Eight supplied inclusive open-charm channels plus constant residual width';
    case '6S'
        info.status='reserved_unavailable'; info.usage='reserved';
        info.mass=4.5777; info.reference_mass=4.612;
        info.blocker=['6S running BW partial widths are unavailable: ', ...
            'DsstarDs1prime(2535) threshold 4.6472 GeV exceeds REF=4.612 GeV ', ...
            'although the source lists 0.35 MeV; DDstar(2S) and DstarD1prime ', ...
            'lack broad-daughter parameters. Registered M=4.5777 GeV.'];
        info.truncation_note=['Known open/gg/gamma-gamma entries sum arithmetically to ', ...
            '25.68120 MeV; this subtotal is not a complete runnable width model'];
end
end
