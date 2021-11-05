<<<<<<< HEAD
function [ELw, EWw, EWd] = get_LW_j(par, cPar, tT, time, f_egg, TC_ah, f_fish) 
%input: 
% par
% cPar
% % -------------------  for embryo ------------
% pars_UE0 % for computing initial_scaled_reserve
% f_egg % f of mother at time of egg laying
% TC_ah % temp of egg incubation
% % -------------------------------------- 
% E_Hh
% E_Hb
% p_Am
% tT % knots for different rearing temperatures through time
% time % time vector
% pars_T % vector with pars for tempcorr (T_A, T_ref)
% 
 vars_pull(par); 
  vars_pull(cPar); 
f = f_fish;

%% calculate the condition at hatching 
% calculate the initial scaled reserve for a specific f of the mother at time of egg laying
 pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f_egg, pars_UE0); % d.cm^2, initial scaled reserve

U_Hh = E_Hh/p_Am; U_Hb = E_Hb/p_Am; 
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
E_h = aUL(2,2) * p_Am *TC_ah; L_h = aUL(2,3);     

%% solve equations for growth and mat from hatching onwards
ELH_h = [E_h/L_h^3 , L_h , E_Hh]; % initial conditions (at hatching), using reserve density!!
% tT = [0 20 65 100; f_j_JH78 f_j_JH78 f f]'; % different temps throughout rearing

options = odeset('Events', @event_bj); % capture events - birth and metamorphosis

pars_lj =  [g, k, l_T, v_Hb, v_Hj ];
[lj, ~,  lb, info ] = get_lj(pars_lj, f); 

if size(tT,1)>1
    tTC = [tT(:,1), tempcorr(tT(:,2),T_ref, T_A)]; %make vector of temp corrections
else
    tTC = tempcorr(tT,T_ref, T_A); % use tempcorr from dataset
    

[t, ELH, te, ye, ie] = ode45(@get_ELH_j, time, ELH_h, options, f, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tTC, lb, lj); % ELH: {J/cm^3, cm}, with {[E], L, H}


% ELH(end,:) = [];
  EL = ELH(:,2); Lj = ye(2,2); % struct length  and struct length at metam; 
  ELw = [EL(EL<Lj) * del_Me; EL(EL>=Lj) *del_M]; 
  EWw = ELH(:,2).^3 .* (1 + ELH(:,1) * w_E/ mu_E/ d_E); % g, wet weight
  EWd = ELH(:,2).^3 * d_V .* (1 + ELH(:,1) * w_E/ mu_E); % g, dry weight (CHECK!)
%   Ww_b = ye(1,2).^3 .* (1 + ye(1,1) * w_E/ mu_E/ d_E); % g, wet weight at birth
%   tT_j = te(2);  % d, time since birth at metam at f and T
%   Ww_j = ye(2,2).^3 .* (1 + ye(2,1) * w_E/ mu_E/ d_E); % g, wet weight at metam

end

function dELH = get_ELH_j(t, ELH, f, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tTC, lb, lj)
  E = ELH(1); % J/cm^3, reserve density [E]
  L = ELH(2); % cm, structural length
  E_H = ELH(3); % J, maturity (cumulated energy invested in maturation)
  
  if length(tTC)>1 
    TC  = spline1(t, tTC);     % find the correct temperature correction at t 
  end
    
  % % need to include: switch at birth (to f != 0 and V1), switch at metamorph (to isomorph)
  % V1 morph : apply correction to  v, {pAm} (and {F_m} and {p_T}); sM = lj/lb (acceleration factor)
  
  if E_H < E_Hb
     sM = 1
     f = 0;
  elseif E_H >= E_Hb && E_H < E_Hj
      l = L/L_m; sM = l/lb
  else
      sM = lj/lb  
  end
  %v = v * sM; p_Am = sM * p_Am; 
  
  vT = sM * v * TC; pT_Am = sM * p_Am * TC; kT_J = k_J *TC; 
  
  dE = (f *  pT_Am - E * vT)/ L;            % J/d.cm^3, change in reserve density d/dt [E]
  e  = E/ E_m;                            % -, scaled reserve density
  rT  = vT * (e/ L - 1/ (L_m*sM))/ (e + g);      % 1/d, specific growth rate
  dL = L * rT/ 3;                         % cm/d, change in structural length d/dt L
  pC  = E * L^3 * (vT/ L - rT);
  dE_H   = max(0, ((1 - kap) * pC - kT_J * E_H) );     % J, change in cumulated energy invested in maturation

  dELH = [dE; dL; dE_H]; % catenate for output
end

function [value,isterminal,direction] = event_bj(t, ELH, f, v, g, ~, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tT, pars_T, lb, lj)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis

  value = [E_Hb; E_Hj] - ELH(3);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end
%


=======
function [ELw, EWw, EWd] = get_LW_j(par, cPar, tT, time, f_egg, TC_ah, f_fish) 
%input: 
% par
% cPar
% % -------------------  for embryo ------------
% pars_UE0 % for computing initial_scaled_reserve
% f_egg % f of mother at time of egg laying
% TC_ah % temp of egg incubation
% % -------------------------------------- 
% E_Hh
% E_Hb
% p_Am
% tT % knots for different rearing temperatures through time
% time % time vector
% pars_T % vector with pars for tempcorr (T_A, T_ref)
% 
 vars_pull(par); 
  vars_pull(cPar); 
f = f_fish;

%% calculate the condition at hatching 
% calculate the initial scaled reserve for a specific f of the mother at time of egg laying
 pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f_egg, pars_UE0); % d.cm^2, initial scaled reserve

U_Hh = E_Hh/p_Am; U_Hb = E_Hb/p_Am; 
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
E_h = aUL(2,2) * p_Am *TC_ah; L_h = aUL(2,3);     

%% solve equations for growth and mat from hatching onwards
ELH_h = [E_h/L_h^3 , L_h , E_Hh]; % initial conditions (at hatching), using reserve density!!
% tT = [0 20 65 100; f_j_JH78 f_j_JH78 f f]'; % different temps throughout rearing

options = odeset('Events', @event_bj); % capture events - birth and metamorphosis

pars_lj =  [g, k, l_T, v_Hb, v_Hj ];
[lj, ~,  lb, info ] = get_lj(pars_lj, f); 

if size(tT,1)>1
    tTC = [tT(:,1), tempcorr(tT(:,2),T_ref, T_A)]; %make vector of temp corrections
else
    tTC = tempcorr(tT,T_ref, T_A); % use tempcorr from dataset
    

[t, ELH, te, ye, ie] = ode45(@get_ELH_j, time, ELH_h, options, f, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tTC, lb, lj); % ELH: {J/cm^3, cm}, with {[E], L, H}


% ELH(end,:) = [];
  EL = ELH(:,2); Lj = ye(2,2); % struct length  and struct length at metam; 
  ELw = [EL(EL<Lj) * del_Me; EL(EL>=Lj) *del_M]; 
  EWw = ELH(:,2).^3 .* (1 + ELH(:,1) * w_E/ mu_E/ d_E); % g, wet weight
  EWd = ELH(:,2).^3 * d_V .* (1 + ELH(:,1) * w_E/ mu_E); % g, dry weight (CHECK!)
%   Ww_b = ye(1,2).^3 .* (1 + ye(1,1) * w_E/ mu_E/ d_E); % g, wet weight at birth
%   tT_j = te(2);  % d, time since birth at metam at f and T
%   Ww_j = ye(2,2).^3 .* (1 + ye(2,1) * w_E/ mu_E/ d_E); % g, wet weight at metam

end

function dELH = get_ELH_j(t, ELH, f, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tTC, lb, lj)
  E = ELH(1); % J/cm^3, reserve density [E]
  L = ELH(2); % cm, structural length
  E_H = ELH(3); % J, maturity (cumulated energy invested in maturation)
  
  if length(tTC)>1 
    TC  = spline1(t, tTC);     % find the correct temperature correction at t 
  end
    
  % % need to include: switch at birth (to f != 0 and V1), switch at metamorph (to isomorph)
  % V1 morph : apply correction to  v, {pAm} (and {F_m} and {p_T}); sM = lj/lb (acceleration factor)
  
  if E_H < E_Hb
     sM = 1
     f = 0;
  elseif E_H >= E_Hb && E_H < E_Hj
      l = L/L_m; sM = l/lb
  else
      sM = lj/lb  
  end
  %v = v * sM; p_Am = sM * p_Am; 
  
  vT = sM * v * TC; pT_Am = sM * p_Am * TC; kT_J = k_J *TC; 
  
  dE = (f *  pT_Am - E * vT)/ L;            % J/d.cm^3, change in reserve density d/dt [E]
  e  = E/ E_m;                            % -, scaled reserve density
  rT  = vT * (e/ L - 1/ (L_m*sM))/ (e + g);      % 1/d, specific growth rate
  dL = L * rT/ 3;                         % cm/d, change in structural length d/dt L
  pC  = E * L^3 * (vT/ L - rT);
  dE_H   = max(0, ((1 - kap) * pC - kT_J * E_H) );     % J, change in cumulated energy invested in maturation

  dELH = [dE; dL; dE_H]; % catenate for output
end

function [value,isterminal,direction] = event_bj(t, ELH, f, v, g, ~, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tT, pars_T, lb, lj)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis

  value = [E_Hb; E_Hj] - ELH(3);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end
%


>>>>>>> b0b68b681feeea3a0aa228ebca2f7ba83cd4b309
