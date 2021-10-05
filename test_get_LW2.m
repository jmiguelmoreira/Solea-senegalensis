function  test_get_LW2() 
% testing
% to do: run_Solea_sen... with 'keyboard' after calculating temp corrections; then evaluate the following lines 

load pars_Solea.mat %loads parameter and compound parameters for Solea_senegalensis, also data and auxData
vars_pull(par); vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
 
TC_ah = tempcorr(C2K(20), T_ref, T_A);
TC_ab = tempcorr(C2K(20), T_ref, T_A);
TC_aj = tempcorr(C2K(20), T_ref, T_A);
 TC_tL2 = tempcorr(C2K(20), T_ref, T_A);


% assume eggs from aquaculture, Manchado data f = f_Man
 % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f_Man, pars_UE0); % d.cm^2, initial scaled reserve
  
   E_0 = U_E0 * p_Am ;          % J, energy in egg -- does this need to be TEMP CORRECTED? 
      
 % HATCH  & BIRTH with dget_aul (DEBtool routine)
 [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h = aUL(2,1)/ TC_ah;                   % d, age at hatch at f and T
  aT_b11 = aUL(3,1)/ TC_ab;                   % d, age at birth at f and T
  
  E_h = aUL(2,2) * p_Am;             % J, energy in reserves at hatch -- TEMP CORRECTED? no
  E_b = aUL(3,2) * p_Am;             % J, energy in reserves at birth -- TEMP CORRECTED? no
 
  L_h = aUL(2,3); Lw_h = L_h/del_Me;        % cm, (strucural) length at hatch
  L_b11 = aUL(3,3); Lw_b11 = L_b11/del_Me;  % cm, (strucural) length at hatch
  
  e_b = (E_b/L_b11^3)/E_m ;
    
 % BIRTH & METAM with my function + events
   EL_h = [E_h/L_h^3 , L_h]; % initial conditions (at hatching), using reserve density!!
%    tT = temp.Taj; % special time-temp vector for Manchado data
    time = [0 5 6 7 21 100]'; temp = [20 20 20 18 16 16]' ; temp = [20 20 20 20 20 20]' ; 
    tT =  [time , temp]; % simple version with constant temp ili changing temp 
    tTC = [tT(:,1), tempcorr(C2K(tT(:,2)),T_ref, T_A)]; %make vector of temp corrections

  options = odeset('AbsTol',1e-9, 'RelTol',1e-9, 'Events',@event_bj);
  pars_lj =  [g, k, l_T, v_Hb, v_Hj ];
  [lj, ~,  lb, info ] = get_lj(pars_lj, f_Man);  % get Lb and Lj for the specified f
  [t, EL, te, ye, ie] = ode45(@get_EL_j, linspace(time(1),time(end),1e3), EL_h, options, f_Man, v, g, E_m, L_m, p_Am, kap, k_J, lb*L_m, lj*L_m, tTC); % ELH: {J/cm^3, cm}, with {[E], L, H}


  aT_b12 = te(1) + aT_h;
%    tT_j1 = te(2) + aT_h - aT_b12; % time since birth at metam
   aT_j1 = te(2) + aT_h ; % age at metam
  
  L_b12 = ye(1,2);
  L_j1 = ye(2,2);
  
  
%% standard with get_tj
 % life cycle 
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_Man);
aT_b2 = tau_b/ k_M/ TC_ab;           % d, age at birth at f and T
tT_j2 = tau_j/ k_M/ TC_aj - aT_b2;
aT_j2 = tau_j/ k_M/ TC_aj ;
rT_B = rho_B * k_M * TC_aj;  % 1/d, von Bert growth rate 

L_b2 = l_b *L_m; 
L_j2 = l_j * L_m;

s_M2 = l_j/l_b ;

%% ----------------- UNI-VARIATE -----------
%t-L age length juveniles
tL2 = [ ...
3.6   0.00001
10  0.3  
20  0.7
36	1.3
45	1.5
59	1.7
66	1.7
75	1.9
82	2.2 
100 4
150 6
300 8
500 10];  % cm, total length at f and T
tL2(:,1) = tL2(:,1) - 3.6;  % correct to time since birth
% units.tL2   = {'d', 'cm'}; label.tL2 = {'time since birth', 'total length'};  % label.tL2 = {'time since hatching', 'total length'};  
% temp.tL2    = C2K(20);  units.temp.tL2 = 'K'; label.temp.tL2 = 'temperature';
% bibkey.tL2 = 'RibeEngr2017';
 

f_RibeEngr = 1;

  % tL2 RibeEngr2017 - version with ODE
  pars_lj =  [g, k, l_T, v_Hb, v_Hj ];
 [lj, ~,  lb, info ] = get_lj(pars_lj, f_RibeEngr);  % get Lb and Lj for the specified f
  EL_b = [f_RibeEngr*E_m , lb*L_m]; % initial conditions (at birth), using reserve density via maternal effect
   options = odeset('AbsTol',1e-15, 'RelTol',1e-10, 'Events',@event_bj); % increase integration sensitivity; capture events - birth and metamorphosis
 [t, EL, te, ye, ie] = ode23(@get_EL_j, [ tL2(:,1)], EL_b, options, f_RibeEngr, v, g, E_m, L_m, p_Am, kap, k_J, lb*L_m, lj*L_m, TC_tL2); % ELH: {J/cm^3, cm}, with {[E], L, H}
%   EL(1,:) = [];
  ELw2 = EL(:,2) / del_M;
  
% tL2 RibeEngr2017 - with get_tj and vonB
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_RibeEngr);
  kT_M2 = TC_tL2 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis
  L_b = (l_b*L_m);
  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj2 = L_b * exp(tL2((tL2(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji2 = L_i - (L_i - L_j) * exp( - rT_B * (tL2((tL2(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  EL2 = [EL_bj2; EL_ji2]; ELw2a = [EL_bj2/del_M; EL_ji2/del_M]; %TL2 is for juveniles all already metamorphosed

  

 
  
  
  

%% output results
fprintf(1, 'Zero variate data: \n' )
fprintf(1, '    ODEs   get_tj  \n' )
fprintf(1, 'aT_h    %2.4f d, --  \n', aT_h)
fprintf(1, 'aT_b    %2.4f d, %2.4f d, %2.4f d \n', aT_b11, aT_b12, aT_b2)
fprintf(1, 'aT_j    %2.4f d, %2.4f d \n', aT_j1, aT_j2)
fprintf(1, 'L_h    %2.4f cm, --  \n', L_h)
fprintf(1, 'L_b    %2.4f cm, %2.4f cm, %2.4f cm \n', L_b11, L_b12, L_b2)
fprintf(1, 'L_j    %2.4f cm, %2.4f cm \n', L_j1, L_j2)
fprintf(1, 'e_b with ODE is %2.4f , and f_Man is %2.4f  \n', e_b, f_Man)
fprintf(1, 'Uni-variate tL data: \n')
i = EL(1,2) == EL2(1)
 [EL(:,2) , EL2]
figure
hold on
plot(tL2(:,1), tL2(:,2), 'o')
plot(tL2(:,1), ELw2, 'k-')
plot(tL2(:,1), ELw2a, 'k--')
legend({'data','ODE', 'vonB'})
% figure
% plot(tL2(:,1), EL(:,1)/E_m, '-')





end

%% subfunctions for solving ODE, important for larva data
function dEL  = get_EL_j(t, EL, f, v, g, E_m, L_m, p_Am, kap, k_J, L_b, L_j, tTC)
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
%   pC  = E * L^3 * (vT/ L - rT);
%   dE_H   = max(0, ((1 - kap) * pC - kT_J * E_H) );     % J, change in cumulated energy invested in maturation

  dEL = [dE; dL]; % catenate for output
end

function [value,isterminal,direction] = event_bj(t, EL, f, v, g, ~, L_m, p_Am, kap, k_J, L_b, L_j, tT)
  % ELH: 3-vector with state variables [E], L, E_H
  % function to find events at hatching and metamorphosis
  
  value = [L_b; L_j] - EL(2);
  isterminal = [0; 0]; % NO stop at life events
  direction = [0; 0];  
end