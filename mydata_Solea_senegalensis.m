function [data, auxData, metaData, txtData, weights] = mydata_Solea_senegalensis

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Pleuronectiformes'; 
metaData.family     = 'Soleidae';
metaData.species    = 'Solea_senegalensis'; 
metaData.species_en = 'Senegalese sole'; 
metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MAE'};
metaData.ecoCode.habitat = {'0jMp', 'jiMcd'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {'Mo'};
metaData.ecoCode.food    = {'bjPz', 'jiCi'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};
metaData.T_typical  = C2K(10); % K, body temp
metaData.data_0     = {'ah'; 'ab'; 'aj'; 'ap'; 'Lh'; 'Lb'; 'Lj'; 'Li'; 'Wwb'; 'Wwp'; 'Wwi'}; 
metaData.data_1     = {'t-L'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Lorna Teal'};    
metaData.date_subm = [2011 08 24];              
metaData.email    = {'lorna.teal@wur.nl'};            
metaData.address  = {'IMARES, IJmuiden'}; 

metaData.curator   = {'Starrlight Augustine'};
metaData.email_cur = {'starrlight@akvaplan.niva.no'}; 
metaData.date_acc  = [2015 08 28]; 



%% set data
% zero-variate data

% data.ah = 2;  units.ah = 'd'; label.ah = 'age at hatch'; bibkey.ah = 'Dinis&Reis1995';  
%   temp.ah = C2K(17);  units.temp.ah = 'K'; label.temp.ah = 'temperature';
data.ah = 1.58; units.ah = 'd'; label.ah = 'age at hatching'; bibkey.ah = 'YufeParr1999'; 
 temp.ah = C2K(19.5) ; units.temp.ah = 'K'; label.temp.ah = 'temperature';
 comment.ah = '~38h from fecundation to hatching rounded to 2';
data.ab = 2+2;      units.ab = 'd';    label.ab = 'age at birth';            bibkey.ab = 'Ribeiro1999'; 
  temp.ab = C2K(17.25); units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.aj = 19;      units.aj = 'd';    label.aj = 'time since fertilization at end of metamorphosis'; bibkey.aj = 'YufeParr1999'; 
   temp.aj = C2K(19.5);  units.temp.aj = 'K'; label.temp.aj = 'temperature';
% data.aj = 13.81+2;     units.aj = 'd';    label.aj = 'age at metamorphosis';    bibkey.aj = 'Fernandez-Diaz2001';
%   temp.aj = C2K(20); units.temp.aj = 'K'; label.temp.aj = 'temperature';
%   comment.aj = 'end of metam for 50% of sample population at f=1 (high live prey density food regime)';
data.ap = 4 * 365;    units.ap = 'd';    label.ap = 'age at puberty'; bibkey.ap = 'Vina2007';
  temp.ap = C2K(17.5);  units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = 'ctemp is the average between 15 and 20 which is the temperature range for the spawning period i.e. the birthday temperature';
  % data.ap = 110;     units.ap = 'd';    label.ap = 'age at puberty';    bibkey.ap = 'Blanco-Vives2011';
%   temp.ap = C2K(20.7); units.temp.ap = 'K'; label.temp.ap = 'temperature';
%   comment.ap = 'constant T, 80% of population mature'; 
data.am = 11 * 365;    units.am = 'd';    label.am = 'life span'; bibkey.am = 'TeixCabr2010';   %-->they cite a thesis by Andrade 1990 
  temp.am = C2K(17.5);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
  comment.am ='temp bibkey.temp.am=Vina2007 page 21 = average yearly temperature of tagus estuary';
% data.am = 26 * 365; units.am = 'd';   label.am = 'life span';                 bibkey.am = 'Deni1990';   
%   temp.am = C2K(10);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
%   comment.am = 'Northsea estimate of Solea solea'; 

% data.Lh  = 0.266;  units.Lh  = 'cm';  label.Lh  = 'total length at hatch';  bibkey.Lh  = 'Dinis&Ribeiro1995'; 
data.Lh = 0.1382;   units.Lh  = 'cm';  label.Lh  = 'total length at hatching';          bibkey.Lh  = 'YufeParr1999'; 
data.Lb  = 0.28;   units.Lb  = 'cm';  label.Lb  = 'total length at birth';  bibkey.Lb  = 'Ribeiro1999';
data.Lj = 0.9176 ; units.Lj = 'cm'; label.Lj = 'total length at end of metamorphosis'; bibkey.Lj  = 'YufeParr1999';  
    comment.Lj = 'length at day 20 post hatching when metamorphosis ends';
% data.Lj  = 0.867;   units.Lj  = 'cm';  label.Lj  = 'total length at metam';  bibkey.Lj  = 'Fernandez-Diaz2001';
%   comment.Lj = 'end of metam for 50% of sample population at f=1 (high live prey density food regime)';
data.Lp  = 30;     units.Lp  = 'cm';  label.Lp  = 'total length at puberty';  bibkey.Lp  = 'fishbase';
data.Li  = 52;   units.Li  = 'cm';  label.Li  = 'ultimate total length';          bibkey.Li  = 'TeixCabr2010';   %-->they cite a thesis by Andrade 1990 
    comment.Li = 'average size between 40 cm';
    %maximal lenght ever recorded is 70 cm, ref FAO
    %maximal lenght ever recorded is 60 (male) cm, ref Fishbase
% data.Li  = 60;    units.Li  = 'cm';  label.Li  = 'ultimate total length';  bibkey.Li  = 'fishbase';

data.Wdh = 33.19e-6;   units.Wdh = 'g';   label.Wdh = 'dry weight at hatching';   bibkey.Wdh = 'YufeParr1999';
data.Wdb = 35.4e-6;   units.Wdb = 'g';   label.Wdb = 'dry weight at birth';   bibkey.Wdb = 'YufeParr1999';
    comment.Wdb = 'mean dry weigth at day 3 ph';
data.Wdj = 1.281e-3 ; units.Wdj = 'ug'; label.Wdj = 'dry weight at END of metamorphosis'; bibkey.Wdj = 'YufeParr1999'; 
    comment.Wdj = 'mean dry weigth at day 17 ph';
% data.Wdj = 9.3e-07;  units.Wdj = 'g';  label.Wdj = 'dry weight at metam';  bibkey.Wdj = 'Fernandez-Diaz2001';
%   comment.Wdj = 'end of metam for 50% of sample population at f=1 (high live prey density food regime)';
data.Wwi = 3000*(60^3/70^3);  units.Wwi = 'g';    label.Wwi = 'ultimate wet weight';      bibkey.Wwi = 'guess';
  comment.Wwi = 'based on ultimate total length difference with S. solea on fishbase';


data.Ri  = 4160;   units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate per day';     bibkey.Ri_35  = 'DiniRibe1999';   
temp.Ri = C2K(18.5);    units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
    comment.Ri = 'Total weight of eggs daily collection of S. senegalensis during the spawning seasons of 1996 and 1997 divided by Wd of egg and 365';
% data.R45  = 1.1e6/365; units.R45  = '#/d'; label.R45  = 'reprod rate at 45 cm'; bibkey.R45  = 'WittGree1995';   
% temp.R45 = C2K(10); units.temp.R45 = 'K'; label.temp.R45 = 'temperature';
%   comment.R45 = 'estimate for Solea solea';
 
% uni-variate data
% t-L data from otolith back-calculation data for SOLEA SOLEA females!(IMARES data)
% data.tLTeal = [ ... % time since birth (a), length (cm)
% 1	8.909404598
% 2	19.91239656
% 3	26.95918876
% 4	31.44085313
% 5	34.39816096
% 6	36.2481578
% 7	37.61539273
% 8	38.32918647
% 9	39.16698018
% 10	39.22343839
% 11	39.87342987
% 12	40.96987132
% 13	42.20480305
% 14	42.43587572
% 15	44.03172643];
% data.tLTeal(:,1) = 365 * data.tLTeal(:,1); % covert a to d
% units.tLTeal   = {'d', 'cm'};  label.tLTeal = {'time since birth', 'total length'};  
% temp.tLTeal    = C2K(10);  units.temp.tLTeal = 'K'; label.temp.tLTeal = 'temperature';
% bibkey.tLTeal = 'Teal2011';

%t-L adult females age from otolith back-calculation data from the field
data.tL3 = [ ... % years, lenght mm
2	215.57
2	220.996
2	224.403
2	255.512
2	259.594
2	228.052
3	223.352
3	231.236
3	277.781
3	303.12
3	275.47
3	309.202
3	234.961
3	294.027
3	305.247
3	298.191
3	244.469
4	291.186
4	307.807
4	295.963
4	325.32
4	302.172
4	286.569
4	330.727
4	318.952
4	266.529
4	282.57
4	243.043
4	271.146
4	342.132
4	247.978
4	313.061
4	353.99
4	252.882
4	349.385
4	275.903
5	332.045
5	337.139
5	376.298
5	365.682
5	384.099
5	412.1
5	389.671
5	233.299
5	403.203
5	283.958
5	296.529
5	303.388
6	339.126
6	332.21
6	410.157
6	400.166
6	343.38
6	350.831
6	300.783
6	359.963
6	305.515
6	405.05
6	354.735
7	406.292
7	471.795
7	370.432
7	401.834
7	445.896
7	303.286
7	320.773
7	348.921
7	354.595
7	310.743
7	433.318
7	315.201
8	450.207
8	454.824
8	362.205
8	400.809
8	418.111];
data.tL3(:,1) = data.tL3(:,1)*365 - data.ab;    %convert to days
data.tL3(:,2) = data.tL3(:,2)/10;               %convert to cm 
units.tL3   = {'d', 'cm'}; label.tL3 = {'time since birth', 'total length'};  % label.tL2 = {'time since hatching', 'total length'};  
temp.tL3    = C2K(17.5);  units.temp.tL3 = 'K'; label.temp.tL3 = 'temperature';
comment.tL3 = 'Animals sampled from Portuguese coast, temperature is a yearly guesstimate'; 
bibkey.tL3 = 'TeixCabr2010';

%tL adults males
data.tL4 = [ ... % age from otolith back-calculation data years, mm data from the field
2	233.962
2	251.589
2	245.539
2	247.991
2	279.577
2	263.969
3	207.986
3	199.777
4	233.497
4	259.505
5	325.365
5	350.513
5	285.742
5	412.068
5	298.682
5	290.193
5	295.453
5	301.484
5	341.795
5	308.151
6	239.235
6	383.742
6	300.083
6	370.292
6	303.687
6	391.944
6	308.482
7	315.825
7	319.964
7	349.01
7	354.582
7	303.662
7	312.227
8	400.656];
data.tL4(:,1) = data.tL4(:,1)*365 - data.ab;    %convert to days and start at birth
data.tL4(:,2) = data.tL4(:,2)/10;               %convert to cm 
units.tL4   = {'d', 'cm'}; label.tL4 = {'time since birth', 'total length'};  % label.tL2 = {'time since hatching', 'total length'};  
temp.tL4 = C2K(17.5);  units.temp.tL4 = 'K'; label.temp.tL4 = 'temperature';
comment.tL4 = 'Animals sampled from Portuguese coast, temperature is a yearly guesstimate'; 
bibkey.tL4 = 'TeixCabr2010';

data.tLA = [ ... % time (d), total length (cm)
0	13.3
0	12
0	12.5
0	13.4
61	17.5
61	17.8
75	16.9
75	16.8
75	18.4
0	11.5
0	10.7
0	13.3
0	12
61	16.0
61	16.7
75	18.3
75	17.1
75	14.8
0	12
0	12.5
0	10.1
0	10.4
61	14.5
61	15.0
75	17.2
75	16.5];
units.tLA   = {'d', 'cm'};  label.tLA = {'time', 'total length'};  
temp.tLA    = C2K(19);  units.temp.tLA = 'K'; label.temp.tLA = 'temperature';
length0.tLA = mean(data.tLA((data.tLA(:,1)==0),2)); units.length0.tLA = 'cm'; label.length0.tLA = 'initial length';
bibkey.tLA = 'expA';
comment.tLA = 'pH 8.0';

data.tLB = [ ... % time (d), total length (cm)
0	14.7
0	16.5
0	13.5
0	14.3
0	16
61	20.1
61	21.1
75	18.4
75	18.9
75	20.3
0	13.8
0	14
0	13.5
0	13.6
61	17.9
61	18.3
75	17.9
75	16.7
75	18.7
0	12.7
0	11.6
0	13
61	16.4
75	16.1
75	18.4];
units.tLB   = {'d', 'cm'};  label.tLB = {'time', 'total length'};  
temp.tLB    = C2K(19);  units.temp.tLB = 'K'; label.temp.tLB = 'temperature';
length0.tLB = mean(data.tLB((data.tLB(:,1)==0),2)); units.length0.tLB = 'cm'; label.length0.tLB = 'initial length';
bibkey.tLB = 'expB';
comment.tLB = 'pH 7.7';

data.tLC = [ ... % time (d), total length (cm)
0	12
0	13
0	12.5
0	11.5
0	12.7
61	17.3
61	16.7
75	16.2
75	20.2
0	13.8
0	11.5
0	12.4
0	13.3
61	15.7
61	17.0
75	18.3
75	17.9
75	17.8
0	12
0	12.5
0	11.7
0	12
0	10.5
61	18.8
61	18.0
75	17.1
75	18.6];
units.tLC   = {'d', 'cm'};  label.tLC = {'time', 'total length'};  
temp.tLC    = C2K(23);  units.temp.tLC = 'K'; label.temp.tLC = 'temperature';
length0.tLC = mean(data.tLC((data.tLC(:,1)==0),2)); units.length0.tLC = 'cm'; label.length0.tLC = 'initial length';
bibkey.tLC = 'expC';
comment.tLC = 'pH 8.0';

data.tLD = [ ... % time (d), total length (cm)
0	12.3
0	13.7
0	13.6
0	12.8
61	18.5
61	17.6
75	19.5
75	19.0
75	17.2
0	13.2
0	11.2
0	13
0	12
0	13.5
61	18.8
61	16.0
75	18.5
75	17.2
75	19.2
0	15.3
0	12.4
0	12.8
0	12.8
61	20.7
61	17.3
75	18.9];
units.tLD   = {'d', 'cm'};  label.tLD = {'time', 'total length'};  
temp.tLD    = C2K(23);  units.temp.tLD = 'K'; label.temp.tLD = 'temperature';
length0.tLD = mean(data.tLD((data.tLD(:,1)==0),2)); units.length0.tLD = 'cm'; label.length0.tLD = 'initial length';
bibkey.tLD = 'expD';
comment.tLD = 'pH 7.7';

data.tLE = [ ... % time (d), total length (cm)
0	33
1	40.2
2	44
3	46
0	37
1	43.7
2	49
3	50
0	43.2
1	48.1
2	48
3	57
0	40.3
1	46.1
2	50
3	51
0	40.5
1	45.8
2	49
3	51.5
0	35.9
1	39.2
2	45
3	46.5
0	35.1
1	40.8
2	44
3	47
0	36
1	40.2
2	42
3	43
0	37.8
1	42.2
2	49
3	50.5
0	39.4
1	44.6
2	49.5
3	48.5
0	32.3
1	37.3
3	44
0	34.9
1	39.9
2	45
3	45
0	33.7
1	40
2	42
3	47.5
0	37.3
0	35.5
1	42.4
2	48
3	50];
data.tLE(:,1) = data.tLE(:,1) * 365; %convert yr to d
units.tLE   = {'d', 'cm'};  label.tLE = {'time', 'total length'};  
temp.tLE    = C2K(19);  units.temp.tLE = 'K'; label.temp.tLE = 'temperature';
length0.tLE = mean(data.tLE((data.tLE(:,1)==0),2)); units.length0.tLE = 'cm'; label.length0.tLE = 'initial length';
bibkey.tLE = 'Lots';
comment.tLE = 'SsEPPO1.20140323';

%t-Ww data
data.tWwCB = [ ... % time since birth (d), length (cm)
67	3.76
94	6.56
89	4.9
98	5.44
143	7.67
91	4.67
99	4.3
226	11.78
100	7.95
104	8.46
120	7.11
162	7.69
212	10.33];
units.tWwCB   = {'d', 'cm'};  label.tWwCB = {'time since birth', 'total length'};  
temp.tWwCB    = C2K(23);  units.temp.tWwCB = 'K'; label.temp.tWwCB = 'temperature';
bibkey.tWwCB = 'CasteloBranco2003';

data.tWwA = [ ... % time (d), wet weight (g)
0	27.88
0	27.2
0	20.95
0	27.26
0	34.19
38	47.6
38	37.32
38	59
38	48.15
38	70.35
61	95.32
61	71.49
75	57.53
75	87.89
75	74.11
0	27.23
0	19.16
0	19.7
0	26.85
0	22.53
25	37.18
25	37.18
25	31.13
25	36.04
25	26.91
61	59.78
61	60.54
75	69.75
75	71.31
75	52.85
0	15.82
0	20.81
0	18.34
0	13.85
11	23.96
11	20.11
11	27.5
11	18.47
61	65.24
61	47.33
75	53.58
75	46.56];
units.tWwA   = {'d', 'g'};  label.tWwA = {'time', 'wet weight'};   
temp.tWwA    = C2K(19);  units.temp.tWwA = 'K'; label.temp.tWwA = 'temperature';
length0.tWwA = mean(data.tLA((data.tLA(:,1)==0),2)); units.length0.tLA = 'cm'; label.length0.tLA = 'initial length';
bibkey.tWwA = 'expA';
comment.tWwA = 'pH 8.0';

data.tWwB = [ ... % time (d), wet weight (g)
0	33.86
0	35.53
0	30.66
0	28.45
0	43.9
11	36.2
11	49.8
11	42.7
11	37.27
11	44.16
61	94.96
61	103.54
75	59.22
75	93.4
75	88.31
0	30.77
0	29.83
0	32.54
0	26.61
0	31.77
25	51.74
25	50.62
25	50.32
25	47.75
25	41.93
38	56.35
38	54.06
38	56.01
38	56.46
38	46.45
61	78.97
61	82
75	79.04
75	79.83
75	79.64
0	24.98
0	20.38
0	18.44
0	24.05
61	69.88
75	58.33
75	74.44];
units.tWwB   = {'d', 'g'};  label.tWwB = {'time', 'wet weight'};   
temp.tWwB    = C2K(19);  units.temp.tWwB = 'K'; label.temp.tWwB = 'temperature';
length0.tWwB = mean(data.tLB((data.tLB(:,1)==0),2)); units.length0.tLB = 'cm'; label.length0.tLB = 'initial length';
bibkey.tWwB = 'expB';
comment.tWwB = 'pH 7.7';

data.tWwC = [ ... % time (d), wet weight (g)
0	30.69
0	29.37
0	26.78
0	21.02
0	24.1
38	55.08
38	67.8
38	61.3
38	45.04
61	77.27
61	101.44
75	70.28
75	121.53
0	26.52
0	18.09
0	22.54
0	28.97
0	25.42
25	44.83
25	32.14
25	30.78
25	31.54
61	78.72
61	87.9
75	79.62
75	76.43
75	74.94
0	23.56
0	21.13
0	16.94
0	19.77
0	15.81
11	29.3
11	32.55
11	24.25
11	26.7
11	22
61	101.18
61	80.48
75	80.36
75	45.67
75	105.41];
units.tWwC   = {'d', 'g'};  label.tWwC = {'time', 'wet weight'};  
temp.tWwC    = C2K(23);  units.temp.tWwC = 'K'; label.temp.tWwC = 'temperature';
length0.tWwC = mean(data.tLC((data.tLC(:,1)==0),2)); units.length0.tLC = 'cm'; label.length0.tLC = 'initial length';
bibkey.tWwC = 'expC';
comment.tWwC = 'pH 8.0';

data.tWwD = [ ... % time (d), wet weight (g)
0	24.9
0	23.34
0	25.23
0	26.95
0	21.03
61	84.71
61	77.44
75	85.83
75	89.91
75	58.73
75	74.95
0	24.33
0	22.25
0	20.29
0	23.96
0	25.33
25	36.29
25	44.32
25	28.71
25	49.39
25	36.94
38	36.27
38	51.24
38	60.24
38	47.18
38	46.26
61	94.04
61	78.32
75	72.64
75	77.86
75	84.36
0	34.19
0	31.52
0	23.83
0	33.81
11	46.78
11	41.74
11	39.4
11	32.48
61	118.18
61	117.48
75	102.33];
units.tWwD   = {'d', 'g'};  label.tWwD = {'time', 'wet weight'};  
temp.tWwD    = C2K(23);  units.temp.tWwD = 'K'; label.temp.tWwD = 'temperature';
length0.tWwD = mean(data.tLD((data.tLD(:,1)==0),2)); units.length0.tLD = 'cm'; label.length0.tLD = 'initial length';
bibkey.tWwD = 'expD';
comment.tWwD = 'pH 7.7';

data.tWwE = [ ... % time (d), wet weight (g)
0	614
1	1255
2	1500
3	1565
0	872
1	1368
2	1780
3	2020
0	821
1	1491
2	2400
3	2500
0	927
1	1696
2	2100
3	2100
0	924
1	1376
2	1700
3	1800
0	526
1	722
2	1250
3	1274
0	518
1	970
2	1360
3	1253
0	699
1	860
2	1100
3	1035
0	689
1	1128
2	1700
3	1775
0	765
1	1258
2	1600
3	1473
0	550
1	896
3	1612
0	619
1	897
2	1160
3	1112
0	628
1	1017
2	1360
3	1371
0	677
0	628
1	1257
2	2000
3	1925];
data.tWwE(:,1) = data.tWwE(:,1) * 365; %convert yr to d
units.tWwE   = {'d', 'g'};  label.tWwE = {'time', 'wet weight'};  
temp.tWwE    = C2K(19);  units.temp.tWwE = 'K'; label.temp.tWwE = 'temperature';
length0.tWwE = mean(data.tLE((data.tLE(:,1)==0),2)); units.length0.tLE = 'cm'; label.length0.tLE = 'initial length';
bibkey.tWwE = 'Lots';
comment.tWwE = 'SsEPPO1.20140323';

%L-Ww data
data.LWw = [ ... % total length (cm), wet weight (g)
13.3	27.88
12	27.2
12.5	20.95
13.4	27.26
17.5	95.32
17.8	71.49
16.9	57.53
16.8	87.89
18.4	74.11
11.5	27.23
10.7	19.16
13.3	19.7
12	26.85
16.0	59.78
16.7	60.54
18.3	69.75
17.1	71.31
14.8	52.85
12	15.82
12.5	20.81
10.1	18.34
10.4	13.85
14.5	65.24
15.0	47.33
17.2	53.58
16.5	46.56
14.7	33.86
16.5	35.53
13.5	30.66
14.3	28.45
16	43.9
20.1	94.96
21.1	103.54
18.4	59.22
18.9	93.4
20.3	88.31
13.8	30.77
14	29.83
13.5	32.54
13.6	26.61
17.9	78.97
18.3	82
17.9	79.04
16.7	79.83
18.7	79.64
12.7	24.98
11.6	20.38
13	18.44
16.4	69.88
16.1	58.33
18.4	74.44
12	30.69
13	29.37
12.5	26.78
11.5	21.02
12.7	24.1
17.3	77.27
16.7	101.44
16.2	70.28
20.2	121.53
13.8	26.52
11.5	18.09
12.4	22.54
13.3	28.97
15.7	78.72
17.0	87.9
18.3	79.62
17.9	76.43
17.8	74.94
12	23.56
12.5	21.13
11.7	16.94
12	19.77
10.5	15.81
18.8	101.18
18.0	80.48
17.1	80.36
18.6	45.67
12.3	24.9
13.7	23.34
13.6	25.23
12.8	26.95
18.5	84.71
17.6	77.44
19.5	85.83
19.0	89.91
17.2	58.73
13.2	24.33
11.2	22.25
13	20.29
12	23.96
13.5	25.33
18.8	94.04
16.0	78.32
18.5	72.64
17.2	77.86
19.2	84.36
15.3	34.19
12.4	31.52
12.8	23.83
12.8	33.81
20.7	118.18
17.3	117.48
18.9	102.33];
units.LWw = {'cm', 'g'}; label.LWw = {'total length', 'wet weight'};
bibkey.LWw = 'exp2019';
comment.LWw = 'ABCD';

data.LWwE = [ ... % total length (cm), wet weight (g)
33	614
40.2	1255
44	1500
46	1565
37	872
43.7	1368
49	1780
50	2020
50.2	1445
51.7	1437
51	1530
54	1528
43.2	821
48.1	1491
48	2400
57	2500
40.3	927
46.1	1696
50	2100
51	2100
40.5	924
45.8	1376
49	1700
51.5	1800
35.9	526
39.2	722
45	1250
46.5	1274
35.1	518
40.8	970
44	1360
47	1253
36	699
40.2	860
42	1100
43	1035
37.8	689
42.2	1128
49	1700
50.5	1775
39.4	765
44.6	1258
49.5	1600
48.5	1473
32.3	550
37.3	896
44	1612
34.9	619
39.9	897
45	1160
45	1112
33.7	628
40	1017
42	1360
47.5	1371
37.3	677
35.5	628
42.4	1257
48	2000
50	1925];
units.LWwE = {'cm', 'g'}; label.LWwE = {'total length', 'wet weight'};
bibkey.LWwE = 'Lots';
comment.LWwE = 'SsEPPO1.20140323 and Ss3EPPO20120422';



%% set weights for all real data
weights = setweights(data, []);

weights.ap = weights.ap * 0;
weights.am = weights.am * 0;
weights.Wdj = weights.Wdj * 0;

weights.tLA = weights.tLA / 9;
weights.tLB = weights.tLB / 9;
weights.tLC = weights.tLC / 9;
weights.tLD = weights.tLD / 9;
weights.tLE = weights.tLE / 3;
weights.tL3 = weights.tL3 / 2;
weights.tL4 = weights.tL4 / 2;
weights.tWwA = weights.tWwA / 9;
weights.tWwB = weights.tWwB / 9;
weights.tWwC = weights.tWwC / 9;
weights.tWwD = weights.tWwD / 9;
weights.tWwE = weights.tWwE / 3;
weights.LWw = weights.LWw / 9;
weights.LWwE = weights.LWwE / 3;

%Adriana Sardi Weights
% %weights.ah = 5 * weights.ah; 
% %weights.ab = 5 * weights.ab; 
% %weights.Wwj0 = 0 * weights.Wwj0; 
% weights.tWd = 10 * weights.tWd;
% %weights.LWw = 10 * weights.LWw;
% % weights.tL2 = 10 * weights.tL2;
% weights.Ri = 20 * weights.Ri;
% weights.Lh = 5 * weights.Lh;
% weights.Li = 10 * weights.Li;
% weights.Lp_f = 5 * weights.Lp_f;
% weights.Lp_m = 5 * weights.Lp_m;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
auxData.length0 = length0;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Links
metaData.links.id_CoL = 'dfce073112a528ed6b275704bfb20337'; % Cat of Life
metaData.links.id_EoL = '46570284'; % Ency of Life
metaData.links.id_Wiki = 'Solea_senegalensis'; % Wikipedia
metaData.links.id_ADW = 'Solea_senegalensis'; % ADW
metaData.links.id_Taxo = '187842'; % Taxonomicon
metaData.links.id_WoRMS = '127159'; % WoRMS
metaData.links.id_fishbase = 'Solea-senegalensis.html'; % fishbase

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Solea_senegalensis}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'expA'; type = 'experiment'; bib = [ ... 
'author = {FISHBUDGET Project}, ' ... 
'year = {2019}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'expB'; type = 'experiment'; bib = [ ... 
'author = {FISHBUDGET Project}, ' ... 
'year = {2019}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'expC'; type = 'experiment'; bib = [ ... 
'author = {FISHBUDGET Project}, ' ... 
'year = {2019}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'expD'; type = 'experiment'; bib = [ ... 
'author = {FISHBUDGET Project}, ' ... 
'year = {2019}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Lots'; type = 'experiment'; bib = [ ... 
'author = {EPPO}, ' ... 
'year = {2020}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Teal2011'; type = 'Misc'; bib = ...
'note = {IMARES data base frisbe}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'CasteloBranco2003'; type = 'Thesis'; bib = [ ...  
'author = {Castelo Branco, Maria Ana}, ' ...
'year = {2003}, ' ...
'title  = {Estudo da producao de linguado (Solea senegalensis Kaup,1858) em tanques de terra}, ' ...
'publisher = {IPIMAR, Instituto Nacional de Investigacao Agraria e das Pescas}, ' ...
'pages = {Table 2.3. (page 26)}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'TeixCabr2010'; type = 'Article'; bib = [ ... 
'author = {﻿Teixeira, Célia M.; Cabral, Henrique N.}, ' ... 
'year   = {2010}, ' ...
'title  = {Comparative analysis of the diet, growth and reproduction of the soles, \textit{Solea solea} and \textit{Solea senegalensis}, occurring in sympatry along the Portuguese coast}, ' ...
'journal= {Journal of the Marine Biological Association of the United Kingdom}, ' ...
'volume = {90}, ' ...
'number = {5}, '...
'doi    = {﻿10.1017/S0025315410000238}, '...
'pages  = {995-1003}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Vina2007'; type = 'Thesis'; bib = [ ... 
'author = {Catarina Vinagre}, ' ... 
'year   = {2007}, ' ...
'title  = {Ecology of the juveniles of the soles \textit{Solea solea} ({L}innaeus, 1758) and \textit{Solea senegalensis} ({K}aup, 1858), in the {T}agus estuary}, ' ...
'journal= {University of Lisbon}, ' ...
'volume = {}, ' ...
'number = {}, '...
'doi    = {}, '...
'pages  = {1-200}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'YufeParr1999'; type = 'Article'; bib = [ ... 
'author = {Yufera, M. Parra, G. Santiago, R. Carrascosa, M.}, ' ... 
'year   = {1999}, ' ...
'title  = {Growth, carbon, nitrogen and caloric content of \textit{Solea senegalensis} (Pisces: Soleidae) from egg fertilization to metamorphosis}, ' ...
'journal= {Marine Biology}, ' ...
'volume = {134}, ' ...
'number = {1}, '...
'doi    = {10.1007/s002270050523}, '...
'pages  = {43-49}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
%still lacks zero-variate references
%
bibkey = 'fishbase'; type = 'Misc'; bib = ...
'howpublished = {\url{http://www.fishbase.org}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
