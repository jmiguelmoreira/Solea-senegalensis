function [prdData, info] = predict_Solea_senegalensis(par, data, auxData)
  %modified the 10 of February by A. Sardi
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  f_field = f ; % f for field and zero-var data
  
 % filters <-- you need the 'abj' specific filter (with s_M)
 filterChecks =   E_Hh > E_Hb || E_Hh <= 0 || f > 1 || f_tL > 1 || f_RibeEngr > 1 || ... % maturity at hatching has to be between 0 and Ehb
     f_tL < 0.1 || f_RibeEngr < 0.1 || f_exp > 1 || ...
     ~reach_birth(g, k, v_Hb, f_tL) || ... % constraint required for reaching birth with that f
     ~reach_birth(g, k, v_Hb, f_RibeEngr) ; %|| ...  % constraint required for reaching birth with that f
   
  %~reach_birth(g, k, v_Hb, f) || ...
 
  % k * v_Hp >= f^3 || ... % constraint constraint required for reaching puberty with f_field
%      ~reach_birth(g, k, v_Hb, f_TeixCabr) || ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_TeixCabr^3 || ... % constraint constraint required for reaching puberty with f_TeixCabr
%      ~reach_birth(g, k, v_Hb, f_TeixCabr2) || ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_TeixCabr2^3 || ... % constraint constraint required for reaching puberty with f_TeixCabr2
%      ~reach_birth(g, k, v_Hb, f_CanaFern1)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern1^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern1
%      ~reach_birth(g, k, v_Hb, f_CanaFern2)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern2^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern2
%      ~reach_birth(g, k, v_Hb, f_CanaFern3)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern3^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern3
%      ~reach_birth(g, k, v_Hb, f_CanaFern4)|| ... % constraint required for reaching birth with that f
%      k * v_Hp >= f_CanaFern4^3 ||...  % constraint constraint required for reaching puberty with f_CanaFern4
%      k * v_Hp >= f_DiniRibe^3 ||...     % constraint constraint required for reaching puberty with f_DineRibe
%     ~reach_birth(g, k, v_Hb, f_DiniRibe); % constraint required for reaching birth with that f
%       

 
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
 
% 0-var data 'ah';'ab';'tj';'aj';'tp';'am';'Lh';'Lb';'Lp';'Li';'Wwh';'Wwb';'Wdh';'Wdb';'Ri'
  % compute temperature correction factors
  TC_ah = tempcorr(temp.ah, T_ref, T_A);
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_tj = tempcorr(temp.tj, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  % univariate data temp corrections
  TC_tL = tempcorr(temp.tL, T_ref, T_A);
  TC_tL2 = tempcorr(temp.tL2, T_ref, T_A);
  TC_tL_f = tempcorr(temp.tL_f, T_ref, T_A);
  TC_tL_m = tempcorr(temp.tL_m, T_ref, T_A);
  TC_tWd = tempcorr(temp.tWd(1), T_ref, T_A);
  TC_tWd2 = tempcorr(temp.tWd2(1), T_ref, T_A);
  TC_tWd_f1 = tempcorr(temp.tWd_f1(1), T_ref, T_A);
  TC_tWd_f2 = tempcorr(temp.tWd_f2(1), T_ref, T_A);
  TC_tWd_f3 = tempcorr(temp.tWd_f3(1), T_ref, T_A);
  TC_tWd_f4 = tempcorr(temp.tWd_f4(1), T_ref, T_A);
  
   %MARE2019
   % temperatures during the experiments:
   TC_LA = tempcorr(temp.tLA, T_ref, T_A);
   TC_WA = tempcorr(temp.tWwA, T_ref, T_A);
   TC_LB = tempcorr(temp.tLB, T_ref, T_A);
   TC_WB = tempcorr(temp.tWwB, T_ref, T_A);
   % average temperatures between birth and start of the experiments:
   TC_bLA = tempcorr(temp2.tLA, T_ref, T_A);
   TC_bLB = tempcorr(temp2.tLB, T_ref, T_A);
   TC_bWA = tempcorr(temp2.tWwA, T_ref, T_A);
   TC_bWB = tempcorr(temp2.tWwB, T_ref, T_A);
%   
  %% % zero-variate data

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_field);
  
  if info ~= 1 % numerical procedure failed
     fprintf('warning: invalid parameter value combination for get_tj \n')
  end
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  % ^-- U_E0 is underpredicted; this gives an E_0 of U_E0 * p_Am * TC_ah = 2.13e-04

  %EGG
   
  E_0 = U_E0 * p_Am ;          % J, energy in egg
  Wd_0 = E_0 * w_E/ mu_E;      % g, egg dry weight 
  V0 = Wd_0/ d_E;             % cm^3, egg volume 
  Lw_0 = (6 * V0/ pi)^(1/3);  % cm, egg diameter
  
  
  
  % HATCH  
 [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h = aUL(2,1)/ TC_ah;                   % d, age at hatch at f and T
  L_h = aUL(2,3);                           % cm, strucural length at hatch
  Lw_h = L_h/del_Me; 
  E_h = aUL(2,2) * p_Am *TC_ah;             % J, energy in reserves at hatch
  Wd_h = (d_V * L_h^3 + w_E/ mu_E * E_h) *1e6; % ug, dry weight at hatch; 

  % BIRTH
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_Me;                % cm, physical length at birth at f
  aT_b = tau_b/ k_M/ TC_ab;           % d, age at birth at f and T
% 
  
  Wd_b = d_V *L_b^3 * (1 + f * ome) *1e6; % ug, dry weight at birth at f 

  

  % metamorphosis decide if using start or end of metam (change del_M accordingly)
  L_j = L_m * l_j;                  % cm, structural length at metam
  Lw_j = L_j/ del_Me;                % cm, physical length at START of metam at f 
  
  tT_j = (tau_j - tau_b) / k_M/ TC_tj;  % d, time since birth at metam
  aT_j = tau_j/TC_aj/k_M;           % d, age at metamorphosis 
  
  %Ww_j = L_j^3 * (1 + f * ome) * 1e6;     % ug, wet weight at metam 
  Wd_j = L_j^3 * d_V * (1 + f * ome) * 1e6; % ug, dry weight at metam 
  
  
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_field);
  
  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  aT_p = tau_p / k_M/ TC_ap;        %d, age at puberty at f and T
  Lw_p = L_p/ del_M;                % cm, physical length at puberty at f
  Ww_p = L_p^3 *(1 + f * ome);        % g, wet weight at puberty 
 

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate physical length at f
  Ww_i = L_i^3 * (1 + f_field * ome);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap, kap_R, g, k_J, k_M, L_T, v, U_Hb, U_Hj, U_Hp];
  [R_i, UE0, Lb, Lj, Lp, info]  =  reprod_rate_j(L_i, f, pars_R);
  
  RT_i = TC_Ri * reprod_rate_j(Lw_i * del_M, f, pars_R); %--> prediction is 0

% life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f_field, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  
% males  
  p_Amm = z_m * p_M/ kap;
  E_mm     = p_Amm/v;             % J/cm^3, [E_m], reserve capacity 
  m_Emm    = y_E_V * E_mm / E_G;   % mol/mol, reserve capacity
  omem     = m_Emm * w_E * d_V/ d_E/ w_V; % -, \omega, contribution of ash free dry mass of reserve to total ash free dry biomass
  gm       = E_G/ kap/ E_mm ;      % -, energy investment ratio
  L_mm     = v/ k_M/ gm;           % cm, maximum length
% L_mm = L_m; % assume males and females have same L_m (and L_i)because no data on ultimate length or weight for males vs females
%   pars_tjm = [g k l_T v_Hb v_Hj v_Hpm];
pars_tjm = pars_tj; % assume maturity threshold for puberty is the same
  [tau_jm, tau_pm, tau_bm, l_jm, l_pm, l_bm, l_im, rho_jm, rho_Bm] = get_tj(pars_tjm, f_field);
%   tT_p_m = (tau_pm - tau_bm)/ k_M/ TC_ap; % d, time since birth at puberty
  Lw_p_m = (L_mm * l_pm)/ del_M;                % cm, total length at puberty at f
   Ww_p_m = (L_mm * l_pm)^3 *(1 + f * omem);        % g, wet weight at puberty 
 
  
%pack to output
  prdData.ah = aT_h;
  prdData.ab = aT_b;
  prdData.aj = aT_j; prdData.aj2 = aT_j; % age at START and END of metam. 
  % Because metamorphosis is modeled as a 'discrete event' rather than something lasting more days, 
  % the predicted value should fall between these two observed values
  prdData.tj = tT_j;prdData.tj2 = tT_j;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.L0 = Lw_0;
  prdData.Lh = Lw_h;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;prdData.Lj2 = Lw_j;
  prdData.Lp_f = Lw_p; 
  prdData.Lp_m = Lw_p_m; %lenght at puberty for females and males
  prdData.Li = Lw_i;
  prdData.Wd0 = Wd_0;

  % Because metamorphosis is modeled as a 'discrete event' rather than something lasting more days, 
  % the predicted value should fall between these two observed values
  prdData.Wwp_f = Ww_p;
  prdData.Wwp_m = Ww_p_m; %wet weight at puberty for females and males
  prdData.Wwi = Ww_i;
  prdData.Wdh = Wd_h;
  prdData.Wdb = Wd_b;
  prdData.Wdj = Wd_j; prdData.Wdj2 = Wd_j; % dry weight at START and END of metam. 
  % Because metamorphosis is modeled as a 'discrete event' rather than something lasting more days, 
  % the predicted value should fall between these two observed values
  prdData.Ri = RT_i;
  prdData.E0 = E_0;
%   
  %% ------------- uni-variate data----------------
  % PARAMETERS for egg
%   pars_UE0 = [V_Hb; g; k_J; k_M; v];
%   [U_E0, L_b, info] = initial_scaled_reserve(f, pars_UE0);
%   [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
%   
%   % Tah 
%   a_h = aUL(2,1);                             % d, age at birth
%   Eah = a_h*ones(length(Tah(:,1)),1) ./ TC_Tah; % d, age at birth temp corrected   
  
%%
%
% % time-length tL
%   
%   [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL);
  
  %t-L RibeSara1999
  %time-length since hatching concatenated to juveniles (from day 12)
  kT_M = k_M * TC_tL;
  rT_B = rho_B * kT_M;  % 1/d, von Bert growth rate   
  rT_j = rho_j * kT_M;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M; % time since *birth* at metamorphosis


  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  
  EL_bj = L_b * exp(tL((tL(:,1)<= tT_j),1) * rT_j/3); % exponential growth as V1-morph
%   EL_hj = L_h * exp (tL((tL(:,1)<= aT_j),1) * rT_j/3);
  EL_ji = L_i - (L_i - L_j) * exp( - rT_B * (tL((tL(:,1) >= tT_j),1) - tT_j)); % cm, expected length at time
  ELw = [EL_bj/del_Me; EL_ji/del_M]; % catenate lengths
  
  % tL2 RibeEngr2017 only juveniles
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_RibeEngr);
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
  
   %tL_f (t-L3) TeixCabr2010 females
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_TeixCabr);
  kT_M2 = TC_tL_f * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis

  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  EL_bj_f = Lw_b * exp(tL_f((tL_f(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji_f = L_i - (L_i - L_j) * exp( - rT_B * (tL_f((tL_f(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw_f = [EL_bj_f; EL_ji_f]/del_M; 
  
   %tL_m (t-L4) TeixCabr2010 males
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tjm, f_TeixCabr);
  kT_M2 = TC_tL_m * k_M;
  rT_B =  rho_B * kT_M2;  % 1/d, von Bert growth rate   
  rT_j =  rho_j * kT_M2;  % 1/d, exponential growth rate
  tT_j = (tau_j - tau_b)/ kT_M2; % time since birth at metamorphosis

  L_j = l_j * L_mm; 
  L_i = l_i * L_mm;
  EL_bj_m = Lw_b * exp(tL_m((tL_m(:,1)<= tT_j),1)  * rT_j/3); % exponential growth as V1-morph
  EL_ji_m = L_i - (L_i - L_j) * exp( - rT_B * (tL_m((tL_m(:,1) > tT_j),1)- tT_j)); % cm, expected length at time
  ELw_m = [EL_bj_m; EL_ji_m]/del_M; %
  
  %t-L MARE2019 A
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp);
  kT_M = k_M * TC_bLA; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; %kT_M, rT_j and rT_B at the average temperature before the start of experiment
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_start = L_i - (L_i - L_j) * exp( - rT_B * (min(tLA(:,1)) - tT_j)); %length at start of experiment
  kT_M = k_M * TC_LA; rT_B = rho_B * kT_M; %kT_M and rT_B at the temperature of the experiment
  L_starti = L_i - (L_i - L_start) * exp( - rT_B * (tLA(:,1) - tLA(1,1))); % cm, expected length at time
  ELwA = L_starti/ del_M;  % cm, total length
  
  %t-L MARE2019 B
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp);
  kT_M = k_M * TC_bLB; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; %kT_M, rT_j and rT_B at the average temperature before the start of experiment    
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_start = L_i - (L_i - L_j) * exp( - rT_B * (min(tLB(:,1)) - tT_j)); %length at start of experiment
  kT_M = k_M * TC_LB; rT_B = rho_B * kT_M; %kT_M and rT_B at the temperature of the experiment
  L_starti = L_i - (L_i - L_start) * exp( - rT_B * (tLB(:,1) - tLB(1,1))); % cm, expected length at time
  ELwB = L_starti/ del_M;  % cm, total length
  
  
  %% % time-wet weight t-Ww
  
  %t-Ww MARE2019 A
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp);
  kT_M = k_M * TC_bWA; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; %kT_M, rT_j and rT_B at the average temperature before the start of experiment     
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_start = L_i - (L_i - L_j) * exp( - rT_B * (min(tWwA(:,1)) - tT_j)); %length at start of experiment
  kT_M = k_M * TC_WA; rT_B = rho_B * kT_M; %kT_M and rT_B at the temperature of the experiment
  L_starti = L_i - (L_i - L_start) * exp( - rT_B * (tWwA(:,1) - tWwA(1,1))); % cm, expected length at time
  EWwA = L_starti.^3 * (1 + f_exp * ome); % g, wet weight
  
  %t-Ww MARE2019 B
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_exp);
  kT_M = k_M * TC_bWB; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; %kT_M, rT_j and rT_B at the average temperature before the start of experiment   
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_start = L_i - (L_i - L_j) * exp( - rT_B * (min(tWwB(:,1)) - tT_j)); %length at start of experiment
  kT_M = k_M * TC_WB; rT_B = rho_B * kT_M; %kT_M and rT_B at the temperature of the experiment
  L_starti = L_i - (L_i - L_start) * exp( - rT_B * (tWwB(:,1) - tWwB(1,1))); % cm, expected length at time
  EWwB = L_starti.^3 * (1 + f_exp * ome); % g, wet weight   
 
 %% % length-weight
 %lenght wet weight fish at 190, 398 and 790 days of age (manchado data)
ELWw = (LWw(:,1) * del_M).^3 * (1 + f_field * ome); 

ELWw_f = (LWw_f(:,1) * del_M).^3 * (1 + f_field * ome); %for females
ELWw_m = (LWw_m(:,1) * del_M).^3 * (1 + f_field * omem); %for males
 
 
%Length dry weight --> laboratory conditions assume ab libitum use f=1
%LWd (OrtiFune2019)
L1 = LWd(LWd(:,1)<data.Lj,1) * del_Me; % for data before metamorphosis
L2 = LWd(LWd(:,1)>=data.Lj,1) * del_M; % for data after metamorphosis
ELWd1 = [L1; L2].^3 * d_V*(1 + f_YufeParr * ome)*1e6; % ug, wet weight 

% here we assume that wga is the same before and after metamorphosis 

%LWd2 (YufeParr1999)
L3 = LWd2(LWd2(:,1)<data.Lj,1) * del_Me; %before metamorphosis
L4 = LWd2(LWd2(:,1)>data.Lj,1) * del_M; %after metamorphosis
ELWd2 = [L3; L4].^3 * d_V* (1 + f_YufeParr * ome)*1e6; % ug, dry weight 


%LWd3 (RibeEngr2017) --> they are all metamorphosed 
ELWd3 = (LWd3(:,1) * del_M).^3 * d_V* (1 + f_RibeEngr * ome)*1e6; % ug, dry weight 

%% % time-dry weight
% tWd (YufeParr1999)-->laboratory conditions assume ab libitum use f=1
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_YufeParr);
  L_b = l_b * L_m; 
  L_j = l_j * L_m; 
  L_i = l_i * L_m;
  
    
  tT_j1 = (tau_j - tau_b)/(k_M * TC_tWd);    % d, time since birth at metamorphosis corrected at 19 degrees for dry weight data
  rT_j = rho_j * (k_M * TC_tWd);  
  rT_B = rho_B * (k_M * TC_tWd);  
  L_bj = L_b * exp(tWd((tWd(:,1) <= tT_j1),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd((tWd(:,1) > tT_j1),1) - tT_j1)); % cm, expected length at time
  EWd = [L_bj; L_jm].^3 * d_V * (1 + f_YufeParr * ome) * 1e6 ;
  
 % tWd2 (ParrYufe2001) --> laboratory conditions assume ab libitum use f=1
  tT_j1 = (tau_j - tau_b)/(k_M * TC_tWd2);    % d, time since birth at metamorphosis corrected at 19 degrees for dry weight data
  rT_j = rho_j * (k_M * TC_tWd2);  
  rT_B = rho_B * (k_M * TC_tWd2);  
  L_bj = L_b * exp(tWd2((tWd2(:,1) <= tT_j1),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd2((tWd2(:,1) > tT_j1),1) - tT_j1)); % cm, expected length at time
  EWd2 = [L_bj; L_jm].^3 * d_V * (1 + f_YufeParr * ome) * 1e6 ;
  
 %tWd_Feeding regimes (CañaFern1999)
 %1 L100
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern1);
rT_j = rho_j * (k_M * TC_tWd_f1); 
rT_B = rho_B * (k_M * TC_tWd_f1); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f1); 
L_bj = L_b * exp(tWd_f1(tWd_f1(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f1(tWd_f1(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_1 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern1 * ome); % ug, dry weight

% %2 L50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern2);
rT_j = rho_j * (k_M * TC_tWd_f2); 
rT_B = rho_B * (k_M * TC_tWd_f2); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f2); 
L_bj = L_b * exp(tWd_f2(tWd_f2(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f2(tWd_f2(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_2 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern2 * ome); % ug, dry weight
  
%   
% %3 L100I50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern3);
rT_j = rho_j * (k_M * TC_tWd_f3); 
rT_B = rho_B * (k_M * TC_tWd_f3); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f3); 
L_bj = L_b * exp(tWd_f3(tWd_f3(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f3(tWd_f3(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_3 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern3 * ome); % ug, dry weight
  
%4 L50I50
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_CanaFern4);
rT_j = rho_j * (k_M * TC_tWd_f4); 
rT_B = rho_B * (k_M * TC_tWd_f4); 
tT_j = (tau_j - tau_b)/ (k_M * TC_tWd_f4); 
L_bj = L_b * exp(tWd_f4(tWd_f4(:,1) < tT_j,1) * rT_j/ 3); % cm length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWd_f4(tWd_f4(:,1) >= tT_j,1) - tT_j));   % cm, length after V1-morph period
EWd_4 = 1e6 * [L_bj; L_jm].^3  * d_V * (1 + f_CanaFern4 * ome); % ug, dry weight

%% 

  % pack to output
%   prdData.Tah = Eah;
  prdData.tL = ELw;
  prdData.tL2 =ELw2;
  prdData.tL_f =ELw_f;
  prdData.tL_m =ELw_m; 
  prdData.tLA = ELwA;
  prdData.tLB = ELwB;
  prdData.tWwA = EWwA;
  prdData.tWwB = EWwB;
  prdData.LWw = ELWw;
  prdData.LWw_f = ELWw_f;
  prdData.LWw_m = ELWw_m;

% prdData.T_dw = EWd_j; 
% prdData.Ttbj = Etbj;

  prdData.LWd = ELWd1;
  prdData.LWd2 = ELWd2;
  prdData.LWd3 = ELWd3;
  prdData.tWd = EWd;
  prdData.tWd2 = EWd2;
  prdData.tWd_f1 = EWd_1;
  prdData.tWd_f2 = EWd_2;
  prdData.tWd_f3 = EWd_3;
  prdData.tWd_f4 = EWd_4;


  
