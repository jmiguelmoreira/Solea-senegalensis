function dEL  = get_EL_j(t, EL, v, g, E_m, L_m, p_Am, L_b, L_j, f, tTC)
  % if food constant (we can work with L_b and L_j from get_tj )

  E = EL(1); % J/cm^3, reserve density [E]
  L = EL(2); % cm, structural length
  
  
  if length(tTC)>1 
    TC  = spline1(t, tTC);     % find the correct temperature correction at t
  else
    TC = tTC;
  end
    
  % V1 morph : apply correction to  v, {pAm} (and {F_m} and {p_T}); sM = lj/lb (acceleration factor)
  
  if L < L_b
     sM = 1;
     f = 0; 
  elseif L < L_j   
      sM = L/L_b;
  else 
      sM = L_j/L_b;
  end
  
  %v = v * sM; p_Am = sM * p_Am;
  vT = sM * v * TC; pT_Am = sM * p_Am * TC; %kT_J = k_J *TC; 
  
  dE = (f *  pT_Am - E * vT)/ L;            % J/d.cm^3, change in reserve density d/dt [E]
  e  = E/ E_m ;                          % -, scaled reserve density
  rT  = vT * (e/ L - 1/ (L_m*sM))/ (e + g);      % 1/d, specific growth rate
  dL = L * rT/ 3;                         % cm/d, change in structural length d/dt L
  dEL = [dE; dL]; % catenate for output
end

function [value,isterminal,direction] = event_bj(t, EL, v, g, ~, L_m, p_Am, L_b, L_j, f, tT)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis
  
  value = [L_b; L_j] - EL(2);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end