function [prdData, info] = predict_Solea_senegalensis(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
if f_Teal > 1
    prdData = []; info = 0; return
end

if f_CB > 1
    prdData = []; info = 0; return
end

if f_exp > 1
    prdData = []; info = 0; return
end

if f_E > 1
    prdData = []; info = 0; return
end

if s_pH < 0
    prdData = []; info = 0; return
end

if E_Hh >= E_Hb
    prdData = []; info = 0; return
end
        
  % compute temperature correction factors
  TC_ah = tempcorr(temp.ah, T_ref, T_A);
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_R45 = tempcorr(temp.R45, T_ref, T_A);
  TC_Teal = tempcorr(temp.tLTeal, T_ref, T_A);
  TC_CB = tempcorr(temp.tWwCB, T_ref, T_A);
  TC_A = tempcorr(temp.tWwA, T_ref, T_A);
  TC_B = tempcorr(temp.tWwB, T_ref, T_A);
  TC_C = tempcorr(temp.tWwC, T_ref, T_A);
  TC_D = tempcorr(temp.tWwD, T_ref, T_A);
  TC_E = tempcorr(temp.tWwE, T_ref, T_A);
  

  % zero-variate data

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
     % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  
    % hatch   
  [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h = aUL(2,1)/ TC_ah;          % d, age at hatch at f and T
  M_Eh = J_E_Am * aUL(2,2);        % mol, reserve at hatch at f
  L_h = aUL(2,3);                  % cm, structural length at f
  Lw_h = L_h/ del_M;               % cm, S-V length at hatch at f
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M;                % cm, total length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f 
  a_b = t_b/ k_M; aT_b = a_b/ TC_ab; % d, age at birth at f and T

  % metam
  L_j = L_m * l_j;                  % cm, structural length at birth at f
  Lw_j = L_j/ del_M;                % cm, total length at birth at f
  Ww_j = L_j^3 * (1 + f * w);
  Wd_j = Ww_j * d_V;
  %Wd_j = L_j^3 * d_V + E_j * w_E/ mu_E; % g, ultimate dry weight (remove d_V for wet weight)
  aT_j = t_j/ k_M/ TC_aj;   % d, age at birth at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  aT_p = t_p/ k_M/ TC_ap;           % d, age at birth at f and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  RT_45 = TC_R45 * reprod_rate_j(45 * del_M, f, pars_R);  % #/d, reproduction rate for 45 cm

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % pack to output
  prdData.ah = aT_h;
  prdData.ab = aT_b;
  prdData.aj = aT_j;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lh = Lw_h;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wdj = Wd_j;
  prdData.Wwi = Ww_i;
  prdData.R45 = RT_45;
  
  % uni-variate data
  
  % time-length for late juveniles and adults at f of 0-var data
  % time-length Teal
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_Teal);
  kT_M = k_M * TC_Teal; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_bj = L_b * exp(tLTeal((tLTeal(:,1) <= tT_j),1) * rT_j/ 3); % cm, struc length
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tLTeal((tLTeal(:,1) > tT_j),1) - tT_j)); % cm, struct length
  ELw = [L_bj; L_jm]/ del_M;  % cm, total length
  
      % time-length A
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M = k_M * TC_A; rT_B = rho_B * kT_M;      
  L_0_p = length0.tLA; L_0 = L_0_p * del_M;
  t_0 = - log((L_i-L_0)/(L_i-L_j)) / rT_B - tLA(1,1); %L_i and rT_B during experiment
  L = L_i - (L_i - L_j) * exp( - rT_B * (tLA(:,1) + t_0)); % cm, struct length
  ELwA = L/ del_M;  % cm, total length
  
      % time-length B
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M_s = k_M * (1+s_pH) * TC_B; rho_B_s = rho_B * (g+f_exp) / (g+f_exp*(1+s_pH)); rT_B_s = rho_B_s * kT_M_s;      
  L_i_s = L_i / (1+s_pH);
  L_0_p = length0.tLB; L_0 = L_0_p * del_M;
  t_0 = - log((L_i_s-L_0)/(L_i_s-L_j)) / rT_B_s - tLB(1,1);
  L = L_i_s - (L_i_s - L_j) * exp( - rT_B_s * (tLB(:,1) + t_0)); % cm, struct length
  ELwB = L/ del_M;  % cm, total length
  
      % time-length C
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M = k_M * TC_C; rT_B = rho_B * kT_M;      
  L_0_p = length0.tLC; L_0 = L_0_p * del_M;
  t_0 = - log((L_i-L_0)/(L_i-L_j)) / rT_B - tLC(1,1);
  L = L_i - (L_i - L_j) * exp( - rT_B * (tLC(:,1) + t_0)); % cm, struct length
  ELwC = L/ del_M;  % cm, total length
  
      % time-length D
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M_s = k_M * (1+s_pH) * TC_D; rho_B_s = rho_B * (g+f_exp) / (g+f_exp*(1+s_pH)); rT_B_s = rho_B_s * kT_M_s;      
  L_i_s = L_i / (1+s_pH);
  L_0_p = length0.tLD; L_0 = L_0_p * del_M;
  t_0 = - log((L_i_s-L_0)/(L_i_s-L_j)) / rT_B_s - tLD(1,1);
  L = L_i_s - (L_i_s - L_j) * exp( - rT_B_s * (tLD(:,1) + t_0)); % cm, struct length
  ELwD = L/ del_M;  % cm, total length
  
  %time wet weight A
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M = k_M * TC_A; rT_B = rho_B * kT_M;      
  L_0_p = length0.tWwA; L_0 = L_0_p * del_M;
  t_0 = - log((L_i-L_0)/(L_i-L_j)) / rT_B - tWwA(1,1);
  L = L_i - (L_i - L_j) * exp( - rT_B * (tWwA(:,1) + t_0)); % cm, struct length
  Ww = L.^3 * (1 + f_exp * w);
  EWwA = Ww; %g, wet weight
  
    %time wet weight B
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M_s = k_M * (1+s_pH) * TC_B; rho_B_s = rho_B * (g+f_exp) / (g+f_exp*(1+s_pH)); rT_B_s = rho_B_s * kT_M_s;      
  L_i_s = L_i / (1+s_pH);
  L_0_p = length0.tWwB; L_0 = L_0_p * del_M;
  t_0 = - log((L_i_s-L_0)/(L_i_s-L_j)) / rT_B_s - tWwB(1,1);
  L = L_i_s - (L_i_s - L_j) * exp( - rT_B_s * (tWwB(:,1) + t_0)); % cm, struct length
  Ww = L.^3 * (1 + f_exp * w);
  EWwB = Ww; %g, wet weight
  
    %time wet weight C
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M = k_M * TC_C; rT_B = rho_B * kT_M;      
  L_0_p = length0.tWwC; L_0 = L_0_p * del_M;
  t_0 = - log((L_i-L_0)/(L_i-L_j)) / rT_B - tWwC(1,1);
  L = L_i - (L_i - L_j) * exp( - rT_B * (tWwC(:,1) + t_0)); % cm, struct length
  Ww = L.^3 * (1 + f_exp * w);
  EWwC = Ww; %g, wet weight
  
    %time wet weight D
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp); %before experiment (no stress factor)
        %before experiment
 %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M_s = k_M * (1+s_pH) * TC_D; rho_B_s = rho_B * (g+f_exp) / (g+f_exp*(1+s_pH)); rT_B_s = rho_B_s * kT_M_s;      
  L_i_s = L_i / (1+s_pH);
  L_0_p = length0.tWwD; L_0 = L_0_p * del_M;
  t_0 = - log((L_i_s-L_0)/(L_i_s-L_j)) / rT_B_s - tWwD(1,1);
  L = L_i_s - (L_i_s - L_j) * exp( - rT_B_s * (tWwD(:,1) + t_0)); % cm, struct length
  Ww = L.^3 * (1 + f_exp * w);
  EWwD = Ww; %g, wet weight
  
      %time length and time wet weight E
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_E); %before experiment (no stress factor)
        %before experiment
  %kT_M = k_M * TC_bPA; rT_B = rho_B * kT_M; %no need
  L_j = L_m * l_j; L_i = L_m * l_i;
        %during experiment
  kT_M = k_M * TC_E; rT_B = rho_B * kT_M;      
  L_0_p = length0.tLE; L_0 = L_0_p * del_M;
  t_0 = - log((L_i-L_0)/(L_i-L_j)) / rT_B - tLE(1,1);
  L = L_i - (L_i - L_j) * exp( - rT_B * (tLE(:,1) + t_0)); % cm, struct length
  Ww = L.^3 * (1 + f_E * w);
  ELwE = L / del_M; %cm, length
  EWwE = Ww; %g, wet weight
  
  %time wet weight Castelo Branco
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CB);
  kT_M = k_M * TC_CB; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_bj = L_b * exp(tWwCB(tWwCB(:,1) < tT_j,1) * rT_j/ 3);
  L_ji = L_i - (L_i - L_j) * exp( - rT_B * (tWwCB(tWwCB(:,1) >= tT_j,1) - tT_j)); % cm, expected length at time
  Ww = [L_bj; L_ji].^3 * (1 + f_CB * w);
  ELwCB = [L_bj; L_ji]/ del_M;  % cm, total length
  EWwCB = Ww; %g, wet weight

      %length-wet weight
  ELWw = (LWw(:,1)*del_M).^3 * (1 + f_exp * w); % g, wet weight
  
    %length-wet weight A
  ELWwE = (LWwE(:,1)*del_M).^3 * (1 + f_E * w); % g, wet weight
  
  % pack to output
  prdData.tLTeal = ELw;
  prdData.tLA = ELwA;
  prdData.tLB = ELwB;
  prdData.tLC = ELwC;
  prdData.tLD = ELwD;
  prdData.tLE = ELwE;
  prdData.tWwA = EWwA;
  prdData.tWwB = EWwB;
  prdData.tWwC = EWwC;
  prdData.tWwD = EWwD;
  prdData.tWwE = EWwE;
  prdData.tWwCB = EWwCB;
  prdData.LWw = ELWw;
  prdData.LWwE = ELWwE;
  
  
end
  