function [ELw, EWw, EWd] = get_LW_j(time, EL_init, pars_ELj, tf, tTC, pars_obs) 
% the idea is to have a function which would give observations as output, 
% calling either the get_ELH_j (if f is not constant and maturity is tracked), 
% or get_EL_j in the simpler case
    
% input: 
% * t : time vector
% * initial conditions: ELH or EL vector
% * pars_LWj : 
% ------ 9-par vector for ELH_j: v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj
% --or-- 7-par vector for EL_j:  v, g, E_m, L_m, p_Am, L_b, L_j
% * tf : a scalar with f value, or a nx2 matrix with time-f pairs
% * tTC : a scalar with TC value, or a nx2 matrix with time-TC pairs
% * pars_obs : parameters for calculating observables: del_Me, del_M, w_E, mu_E, d_E, d_V

v = pars_ELj(1);
g = pars_ELj(2);
E_m = pars_ELj(3);
L_m = pars_ELj(4);
p_Am = pars_ELj(5);

if length(pars_ELj) == 9 % more complicated option
    
    kap = pars_ELj(6);
    k_J = pars_ELj(7);
    E_Hb = pars_ELj(8);
    E_Hj = pars_ELj(9);
    
  options = odeset('AbsTol',1e-9, 'RelTol',1e-9, 'Events',@event_bjH);
     
[t, EL, te, ye, ie] = ode45(@get_ELH_j, time, EL_init, options, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tf, tTC);

else % f is constant
    L_b = pars_ELj(6);
    L_j = pars_ELj(7);
    f = tf; 

    options = odeset('AbsTol',1e-9, 'RelTol',1e-9, 'Events',@event_bj);
    [t, EL, te, ye, ie] = ode45(@get_EL_j, time, EL_init(1:2), options, v, g, E_m, L_m, p_Am, L_b, L_j, f, tTC);

end

% keyboard

%% unpack auxiliary pars
del_Me = pars_obs(1);
del_M = pars_obs(2);
w_E = pars_obs(3);
mu_E = pars_obs(4);
d_E = pars_obs(5);
d_V = pars_obs(6);

% EL(end,:) = []; % 
  L = EL(:,2); Lj = ye(1,2); % struct length  and struct length at metam; 
  ELw = [L(L<Lj) * del_Me; L(L>=Lj) *del_M]; 
  EWw = L.^3 .* (1 + EL(:,1) * w_E/ mu_E/ d_E); % g, wet weight
  EWd = L.^3 * d_V .* (1 + EL(:,1) * w_E/ mu_E); % g, dry weight (CHECK!)
%   Ww_b = ye(1,2).^3 .* (1 + ye(1,1) * w_E/ mu_E/ d_E); % g, wet weight at birth
%   tT_j = te(2);  % d, time since birth at metam at f and T
%   Ww_j = ye(2,2).^3 .* (1 + ye(2,1) * w_E/ mu_E/ d_E); % g, wet weight at metam

end

%

function [value,isterminal,direction] = event_bjH(t, ELH, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, f, tTC)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis
  
  value = [E_Hb; E_Hj] - ELH(3);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end

function [value,isterminal,direction] = event_bj(t, EL, v, g, ~, L_m, p_Am, L_b, L_j, f, tT)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis
  
  value = [L_b; L_j] - EL(2);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end