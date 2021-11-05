%% (sub)functions for solving ODE, important for larva data if change in food and temp!
function dELH  = get_ELH_j(t, ELH, tf, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tTC)
  persistent Lb Lj % keeps values in Matlab memory within this function!

  E = ELH(1); % J/cm^3, reserve density [E]
  L = ELH(2); % cm, structural length
  E_H = ELH(3); % J, maturity (cumulated energy invested in maturation)
  
  
  if length(tTC)>1 
    TC  = spline1(t, tTC);     % find the correct temperature correction at t
  else
      TC = tTC;
  end
    
  if length(tf)>1 
    f  = spline1(t, tf);     % find the correct temperature correction at t
  else
      f = tf;
  end

  % V1 morph : apply correction to  v, {pAm} (and {F_m} and {p_T}); sM = lj/lb (acceleration factor)
  
  if E_H < E_Hb
     sM = 1;
     f = 0; 
     Lb = L; 
  elseif E_H >= E_Hb && E_H < E_Hj   
      sM = L/Lb;
      Lj = L; 
  else % E_H >= E_H
      sM = Lj/Lb;
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

function [value,isterminal,direction] = event_bj(t, ELH, f, v, g, E_m, L_m, p_Am, kap, k_J, E_Hb, E_Hj, tTC, tf)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis
  
  value = [E_Hb; E_Hj] - ELH(3);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end
%