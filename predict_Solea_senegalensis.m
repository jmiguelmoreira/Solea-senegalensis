function [prdData, info] = predict_Solea_senegalensis(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
% if f_Teal > 1
%     prdData = []; info = 0; return
% end

if f_CB > 1
    prdData = []; info = 0; return
end

if f_exp > 1
    prdData = []; info = 0; return
end

if f_E > 1
    prdData = []; info = 0; return
end

% if f_34 > 1
%     prdData = []; info = 0; return
% end

% if s_pH < 0
%     prdData = []; info = 0; return
% end

if E_Hh >= E_Hb
    prdData = []; info = 0; return
end

 filterChecks =   E_Hh > E_Hb || E_Hh <= 0 || f_field > 1 || f_tL > 1 || f_tL2 > 1 || ... % maturity at hatching has to be between 0 and Ehb
     f_tL < 0.1 || f_tL2 < 0.1 || ...
     ~reach_birth(g, k, v_Hb, f_tL) || ... % constraint required for reaching birth with that f
     ~reach_birth(g, k, v_Hb, f_tL2) ; %|| ...  % constraint required for reaching birth with that f
if filterChecks  
    info = 0;
    prdData = {}; return
end 
        
  % compute temperature correction factors
  TC_ah = tempcorr(temp.ah, T_ref, T_A);
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_tj = tempcorr(temp.tj, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  %TC_Teal = tempcorr(temp.tLTeal, T_ref, T_A);
  TC_CB = tempcorr(temp.tLCB, T_ref, T_A);
  TC_A = tempcorr(temp.tWwA, T_ref, T_A);
  TC_B = tempcorr(temp.tWwB, T_ref, T_A);
  TC_C = tempcorr(temp.tWwC, T_ref, T_A);
  TC_D = tempcorr(temp.tWwD, T_ref, T_A);
  TC_E = tempcorr(temp.tWwE, T_ref, T_A);
  
    T_pars=[T_A, T_L, T_H, T_AL, T_AH];
  TC_tL = tempcorr(temp.tL1, T_ref, T_pars);
  TC_tL2 = tempcorr(temp.tL2, T_ref, T_pars);
  TC_tL3 = tempcorr(temp.tL3, T_ref, T_pars);
  TC_tL4 = tempcorr(temp.tL4, T_ref, T_pars);
  TC_tWd = tempcorr(temp.tWd(1), T_ref, T_pars);
  TC_tWd2 = tempcorr(temp.tWd2(1), T_ref, T_pars);
  TC_tWd_f1 = tempcorr(temp.tWd_f1(1), T_ref, T_pars);
  TC_tWd_f2 = tempcorr(temp.tWd_f2(1), T_ref, T_pars);
  TC_tWd_f3 = tempcorr(temp.tWd_f3(1), T_ref, T_pars);
  TC_tWd_f4 = tempcorr(temp.tWd_f4(1), T_ref, T_pars);
  TC_3 = tempcorr(temp.tL3, T_ref, T_A);
  TC_4 = tempcorr(temp.tL4, T_ref, T_A);
  

  % zero-variate data

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
     % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  E_0 = U_E0 * p_Am ;          % J, energy in egg
  Wd_0 = E_0 * w_E/ mu_E;      % g, egg dry weight 
  
    % hatch   
  [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h = aUL(2,1)/ TC_ah;          % d, age at hatch at f and T
  M_Eh = J_E_Am * aUL(2,2);        % mol, reserve at hatch at f
  L_h = aUL(2,3);                  % cm, structural length at f
  Lw_h = L_h/ del_Me;               % cm, S-V length at hatch at f
  Ww_h = L_h^3 * (1 + f * w);
%   Wd_h = Ww_h * d_V;
  E_h = aUL(2,2) * p_Am *TC_ah;             % J, energy in reserves at hatch
  Wd_h = (d_V * L_h^3 + w_E/ mu_E * E_h); % ug, dry weight at hatch; 
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_Me;                % cm, total length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f 
%   Wd_b = Ww_b * d_V;
  Wd_b = d_V *L_b^3 * (1 + f * w); % ug, dry weight at birth at f 
  a_b = t_b/ k_M; aT_b = a_b/ TC_ab; % d, age at birth at f and T

  % metam
  L_j = L_m * l_j;                  % cm, structural length at birth at f
  Lw_j = L_j/ del_Me;                % cm, total length at birth at f
  Ww_j = L_j^3 * (1 + f * w);
%   Wd_j = Ww_j * d_V;
  %Wd_j = L_j^3 * d_V + E_j * w_E/ mu_E; % g, ultimate dry weight (remove d_V for wet weight)
  Wd_j = L_j^3 * d_V * (1 + f * w); % ug, dry weight at metam 
  aT_j = t_j/ k_M/ TC_aj;   % d, age at birth at metam
  tT_j = (t_j - t_b) / k_M/ TC_tj;  % d, time since birth at metam

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
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
  RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);          % #/d, max reproduction rate 

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
  prdData.tj = tT_j;
  prdData.Lh = Lw_h;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;
  prdData.Lp_f = Lw_p;  prdData.Lp_m = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wd0 = Wd_0;
  prdData.Wdh = Wd_h;
  prdData.Wdb = Wd_b;
  prdData.Wdj = Wd_j;
  prdData.Wwp_f = Ww_p;prdData.Wwp_m = Ww_p;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  prdData.E0 = E_0;
  
  % uni-variate data
  
  % time-length for late juveniles and adults at f of 0-var data
  % time-length Teal
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_Teal);
%   kT_M = k_M * TC_Teal; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
%   L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
%   L_bj = L_b * exp(tLTeal((tLTeal(:,1) <= tT_j),1) * rT_j/ 3); % cm, struc length
%   L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tLTeal((tLTeal(:,1) > tT_j),1) - tT_j)); % cm, struct length
%   ELw = [L_bj; L_jm]/ del_M;  % cm, total length
  
%     %time-length 3
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_34);
%   kT_M = k_M * TC_3; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
%   L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
%   L_bj = L_b * exp(tL3((tL3(:,1) <= tT_j),1) * rT_j/ 3); % cm, struc length
%   L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL3((tL3(:,1) > tT_j),1) - tT_j)); % cm, struct length
%   ELw3 = [L_bj; L_jm]/ del_M;  % cm, total length
%   
%       %time-length 4
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_34);
%   kT_M = k_M * TC_4; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
%   L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
%   L_bj = L_b * exp(tL4((tL4(:,1) <= tT_j),1) * rT_j/ 3); % cm, struc length
%   L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL4((tL4(:,1) > tT_j),1) - tT_j)); % cm, struct length
%   ELw4 = [L_bj; L_jm]/ del_M;  % cm, total length

    %time length Castelo Branco
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CB);
  kT_M = k_M * TC_CB; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_bj = L_b * exp(tLCB(tLCB(:,1) < tT_j,1) * rT_j/ 3);
  L_ji = L_i - (L_i - L_j) * exp( - rT_B * (tLCB(tLCB(:,1) >= tT_j,1) - tT_j)); % cm, expected length at time
  Ww = [L_bj; L_ji].^3 * (1 + f_CB * w);
  ELwCB = [L_bj; L_ji]/ del_M;  % cm, total length
  
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

      %length-wet weight
  ELWw = (LWw(:,1)*del_M).^3 * (1 + f_exp * w); % g, wet weight
  
    %length-wet weight A
  ELWwE = (LWwE(:,1)*del_M).^3 * (1 + f_E * w); % g, wet weight
  
  
  %% Adriana %%
  %time-length tL
  
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL);
  
  %t-L RibeSara1999
  %time-length since hatching concatenated to juveniles (from day 12)
  kT_M = k_M * TC_tL;
  rT_B = rho_B * kT_M;  % 1/d, von Bert growth rate   
  rT_j = rho_j * kT_M;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M; % time since *birth* at metamorphosis


  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  
  EL_bj = L_b * exp(tL1((tL1(:,1)<= tT_j),1) * rT_j/3); % exponential growth as V1-morph
%   EL_hj = L_h * exp (tL((tL(:,1)<= aT_j),1) * rT_j/3);
  EL_ji = L_i - (L_i - L_j) * exp( - rT_B * (tL1((tL1(:,1) >= tT_j),1) - tT_j)); % cm, expected length at time
  ELw1 = [EL_bj/del_Me; EL_ji/del_M]; % catenate lengths
  
  % tL2 RibeEngr2017 only juveniles
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL2);
  kT_M2 = TC_tL2 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis
%   aT_j = tau_j/ kT_M2 ; % time since *fertilization* at metamorphosis --> closer to hatching (which is where the dataset starts!)
  
  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj2 = Lw_b * exp(tL2((tL2(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji2 = L_i - (L_i - L_j) * exp( - rT_B * (tL2((tL2(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw2 = [EL_bj2; EL_ji2]/del_M; %TL2 is for juveniles all already metamorphosed
  
   %t-L3 TeixCabr2010 females
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_TeixCabr);
  kT_M2 = TC_tL3 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis

  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj3 = Lw_b * exp(tL3((tL3(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji3 = L_i - (L_i - L_j) * exp( - rT_B * (tL3((tL3(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw3 = [EL_bj3; EL_ji3]/del_M; 
  
   %t-L4 TeixCabr2010 males
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_TeixCabr2);
  kT_M2 = TC_tL4 * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis

  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj4 = Lw_b * exp(tL4((tL4(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji4 = L_i - (L_i - L_j) * exp( - rT_B * (tL4((tL4(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw4 = [EL_bj4; EL_ji4]/del_M; %
  
  % length-weight
 %lenght wet weight fish at 190, 398 and 790 days of age (manchado data)
ELWwMan = (LWwMan(:,1) * del_M).^3 * (1 + f_field * w); 

% ELWw2 = (LWw2(:,1) * del_M).^3 * (1 + f_Man * w); %for females
% ELWw3 = (LWw3(:,1) * del_M).^3 * (1 + f_Man * w); %for males
 
 
%Length dry weight --> laboratory conditions assume ab libitum use f=1
%LWd (OrtiFune2019)
L1 = LWd1(LWd1(:,1)<data.Lj,1) * del_Me; % for data before metamorphosis
L2 = LWd1(LWd1(:,1)>=data.Lj,1) * del_Me; % for data after metamorphosis
ELWd1 = [L1; L2].^3 * d_V*(1 + f * w)*1e6; % ug, wet weight 

% here we assume that wga is the same before and after metamorphosis 

%LWd2 (RibeEngr2017) --> they are all metamorphosed 
%L3 = LWd2(LWd2(:,1)<data.Lj,1) * del_Me; %before metamorphosis
L4 = LWd2(LWd2(:,1)>data.Lj,1) * del_Me; %after metamorphosis
ELWd2 = L4.^3 * d_V* (1 + f * w)*1e6; % ug, dry weight 

%LWd3 (YufeParr1999)
L5 = LWd3(LWd3(:,1)<data.Lj,1) * del_Me; %before metamorphosis
L6 = LWd3(LWd3(:,1)>data.Lj,1) * del_M; %after metamorphosis
ELWd3 = [L5; L6].^3 * d_V* (1 + f * w)*1e6; % ug, dry weight 

% time-dry weight
% tWd (YufeParr1999)-->laboratory conditions assume ab libitum use f=1
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  L_b = l_b * L_m; 
  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  
    
  tT_j1 = (tau_j - tau_b)/(k_M * TC_tWd);    % d, time since birth at metamorphosis corrected at 19 degrees for dry weight data
  rT_j = rho_j * (k_M * TC_tWd);  
  rT_B = rho_B * (k_M * TC_tWd);  
  L_bj = L_b * exp(tWd((tWd(:,1) <= tT_j1),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd((tWd(:,1) > tT_j1),1) - tT_j1)); % cm, expected length at time
  EWd = [L_bj; L_jm].^3 * d_V * (1 + f * w) * 1e6 ;
  
 % tWd2 (ParrYufe2001) --> laboratory conditions assume ab libitum use f=1
  tT_j1 = (tau_j - tau_b)/(k_M * TC_tWd2);    % d, time since birth at metamorphosis corrected at 19 degrees for dry weight data
  rT_j = rho_j * (k_M * TC_tWd2);  
  rT_B = rho_B * (k_M * TC_tWd2);  
  L_bj = L_b * exp(tWd2((tWd2(:,1) <= tT_j1),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd2((tWd2(:,1) > tT_j1),1) - tT_j1)); % cm, expected length at time
  EWd2 = [L_bj; L_jm].^3 * d_V * (1 + f * w) * 1e6 ;
  
 %tWd_Feeding regimes (CañaFern1999)
 %1 L100
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern1);
rT_j = rho_j * (k_M * TC_tWd_f1); 
rT_B = rho_B * (k_M * TC_tWd_f1); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f1); 
L_bj = L_b * exp(tWd_f1(tWd_f1(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f1(tWd_f1(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_1 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern1 * w); % ug, dry weight

% %2 L50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern2);
rT_j = rho_j * (k_M * TC_tWd_f2); 
rT_B = rho_B * (k_M * TC_tWd_f2); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f2); 
L_bj = L_b * exp(tWd_f2(tWd_f2(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f2(tWd_f2(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_2 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern2 * w); % ug, dry weight
  
%   
% %3 L100I50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern3);
rT_j = rho_j * (k_M * TC_tWd_f3); 
rT_B = rho_B * (k_M * TC_tWd_f3); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f3); 
L_bj = L_b * exp(tWd_f3(tWd_f3(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f3(tWd_f3(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_3 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern3 * w); % ug, dry weight
  
%4 L50I50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern4);
rT_j = rho_j * (k_M * TC_tWd_f4); 
rT_B = rho_B * (k_M * TC_tWd_f4); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f4); 
L_bj = L_b * exp(tWd_f4(tWd_f4(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f4(tWd_f4(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_4 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern4 * w); % ug, dry weight
  
  
  
  % pack to output
  %prdData.tLTeal = ELw;
  prdData.tLCB = ELwCB;
  prdData.tLA = ELwA;
  prdData.tLB = ELwB;
  prdData.tLC = ELwC;
  prdData.tLD = ELwD;
  prdData.tLE = ELwE;
%   prdData.tL3 = ELw3;
%   prdData.tL4 = ELw4;
  prdData.tWwA = EWwA;
  prdData.tWwB = EWwB;
  prdData.tWwC = EWwC;
  prdData.tWwD = EWwD;
  prdData.tWwE = EWwE;
  prdData.LWw = ELWw;
  prdData.LWwE = ELWwE;
  
  %Adriana
  prdData.tL1 = ELw1;
  prdData.tL2 =ELw2;
  prdData.tL3 =ELw3;
  prdData.tL4 =ELw4; 
  prdData.LWwMan = ELWwMan;
  prdData.LWd1 = ELWd1;
  prdData.LWd2 = ELWd2;
  prdData.LWd3 = ELWd3;
  prdData.tWd = EWd;
  prdData.tWd2 = EWd2;
  prdData.tWd_f1 = EWd_1;
  prdData.tWd_f2 = EWd_2;
  prdData.tWd_f3 = EWd_3;
  prdData.tWd_f4 = EWd_4;
  
  
end
  