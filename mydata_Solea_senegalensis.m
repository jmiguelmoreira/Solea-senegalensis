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

data.ah = 2;  units.ah = 'd'; label.ah = 'age at hatch'; bibkey.ah = 'Dinis&Reis1995';  
  temp.ah = C2K(17);  units.temp.ah = 'K'; label.temp.ah = 'temperature';
data.ab = 2+2;      units.ab = 'd';    label.ab = 'age at birth';            bibkey.ab = 'Ribeiro1999'; 
  temp.ab = C2K(17.25); units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.aj = 9.4+2;     units.aj = 'd';    label.aj = 'age at metamorphosis';    bibkey.aj = 'Fernandez-Diaz2001';
  temp.aj = C2K(20); units.temp.aj = 'K'; label.temp.aj = 'temperature';
  comment.aj = 'start of metam for 50% of sample population at f=1 (high live prey density food regime)';
data.ap = 110;     units.ap = 'd';    label.ap = 'age at puberty';    bibkey.ap = 'Blanco-Vives2011';
  temp.ap = C2K(20.7); units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = 'constant T, 80% of population mature'; 
data.am = 26 * 365; units.am = 'd';   label.am = 'life span';                 bibkey.am = 'Deni1990';   
  temp.am = C2K(10);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
  comment.am = 'Northsea estimate of Solea solea'; 

data.Lh  = 0.266;  units.Lh  = 'cm';  label.Lh  = 'total length at hatch';  bibkey.Lh  = 'Dinis&Ribeiro1995'; 
data.Lb  = 0.28;   units.Lb  = 'cm';  label.Lb  = 'total length at birth';  bibkey.Lb  = 'Ribeiro1999';
data.Lj  = 0.59;   units.Lj  = 'cm';  label.Lj  = 'total length at metam';  bibkey.Lj  = 'Fernandez-Diaz2001';
  comment.Lj = 'start of metam for 50% of sample population at f=1 (high live prey density food regime)';
data.Lp  = 30;     units.Lp  = 'cm';  label.Lp  = 'total length at puberty';  bibkey.Lp  = 'fishbase';
data.Li  = 60;    units.Li  = 'cm';  label.Li  = 'ultimate total length';  bibkey.Li  = 'fishbase';

%data.Wdj = 2.8e-07;  units.Wdj = 'g';  label.Wdj = 'dry weight at metam';  bibkey.Wdj = 'Fernandez-Diaz2001';
  %comment.Wdj = 'start of metam for 50% of sample population at f=1 (high live prey density food regime)';
data.Wwi = 3000*(60^3/70^3);  units.Wwi = 'g';    label.Wwi = 'ultimate wet weight';      bibkey.Wwi = 'guess';
  comment.Wwi = 'based on ultimate total length difference with S. solea on fishbase';

data.R45  = 1.1e6/365; units.R45  = '#/d'; label.R45  = 'reprod rate at 45 cm'; bibkey.R45  = 'WittGree1995';   
temp.R45 = C2K(10); units.temp.R45 = 'K'; label.temp.R45 = 'temperature';
  comment.R45 = 'estimate for Solea solea';
 
% uni-variate data
% t-L data from otolith back-calculation data for females!(IMARES data)
data.tL = [ ... % time since birth (a), length (cm)
1	8.909404598
2	19.91239656
3	26.95918876
4	31.44085313
5	34.39816096
6	36.2481578
7	37.61539273
8	38.32918647
9	39.16698018
10	39.22343839
11	39.87342987
12	40.96987132
13	42.20480305
14	42.43587572
15	44.03172643];
data.tL(:,1) = 365 * data.tL(:,1); % covert a to d
units.tL   = {'d', 'cm'};  label.tL = {'time since birth', 'total length'};  
temp.tL    = C2K(10);  units.temp.tL = 'K'; label.temp.tL = 'temperature';
bibkey.tL = 'Teal2011';

  
%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
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
%lacks references
%
bibkey = 'fishbase'; type = 'Misc'; bib = ...
'howpublished = {\url{http://www.fishbase.org}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
