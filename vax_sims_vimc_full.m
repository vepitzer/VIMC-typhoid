%% Simulations for the big burden model

% Pull in burden dynamic link, which includes country data and agedist
%load('./burden_dynamic_link_08Dec2017b.mat') %IHMEs is the same
load('./burden_dynamic_link_Yale.mat');
load('model_samples_ALL.mat');
clear R0country R0modelcasescountry repcountry avgage_country

% Import list of VIMC countries
country_list = readtable('typhoid-93-list.csv','ReadVariableNames',true);

vimc_country_names=country_list.Country; %gavi73.country_name; %{'Congo, the Democratic Republic of the','Ethiopia','India','Nigeria','Pakistan'};  
vimc_countries=country_list.ISO3; %gavi73.ISO3; %{'COD','ETH','IND','NGA','PAK'};  

% Import list of 145 countries
all_countries = readtable('countrycodes_iso.txt','ReadVariableNames',true);

%% Bring in coverage data
% Previous data from Gavi
%cov_unc = readtable('Gavi_cov_unc.csv');
%cov_con = readtable('Gavi_cov_con.csv');

%coverage for VIMC countries
vimc_cov = readtable('coverage_201910gavi-5_typhoid-routine-default.csv','ReadVariableNames',true);

%cov_unc(cov_unc.Year<2019, :) = [];
%cov_con(cov_con.Year<2019, :) = [];
vimc_cov(vimc_cov.year<2000, :) = [];

% Make list of the Gavi73 countries only
gavi73 = table(unique(cov_unc.ISO3));
gavi73.Properties.VariableNames = {'ISO3'};
gavi73.cn = nan(73,1); 

iso_lmic = iso;
iso_lmic(strcmp(iso_lmic.countryiso, 'TWN'),:) = [];
for i=1:73
   gavi73.cn(i) = find(strcmp(gavi73.ISO3(i), iso_lmic.countryiso));
end
for i=1:73
    x = find(strcmp(gavi73.ISO3{i},country_list.ISO3));
    gavi73.country_name(i) = country_list.Country(x);
end

ncountry=length(vimc_countries);

% Bring in population data
input_pop_data = readtable('201910gavi-5_dds-201910_2_int_pop_both.csv','ReadVariableNames',true);
input_pop_data(:,5:121)=[];
vimc_pop_mat=zeros(101,101,ncountry);

% Age-independent fixed parameters
params.delta = 1/4;
params.alpha = 0.01; % 0.01;
params.omega = 1/104; % -log(1-.25)/52;
% params.epsilon = 2.68*10^(-11);

% Age-specific fixed parameters
agepar.u = [1/(0.75*52); 1/(1.25*52); 1/(3*52); 1/520; 1/520; 0];
agepar.theta = [0.003; 0.003; 0.003; 0.003; 0.021; 0.021];
agepar.theta2 = zeros(6,1);

% Birth rates
mub = [36.6; 23.6; 15.0]; 
% ESA estimate of birth rate in low income: 36.6 for 2010-15
% ESA estimate of birth rate in middle income: 19.7 for 2010-15
% ESA estimate of birth rate in lower-middle income: 23.6 for 2010-15
% ESA estimate of birth rate in upper-middle income: 15.0 for 2010-15

% These are only these ages: 0-4, 5-9, 10-14, 15-20, 20-25. The rest are
% 1-sum(distribution)
low_inc = [0.158914104, 0.14116809, 0.124802931, 0.108397066, 0.091790805];
lmid_inc = [0.108066773, 0.103190497, 0.098970907, 0.094718397, 0.090066389];
umid_inc = [0.072224644, 0.068942809, 0.066816143, 0.068359158, 0.080071929];
% mid_inc = [0.091922178	0.087764088	0.084487406	0.08284539	0.085564619];
% Download from ESA. https://esa.un.org/unpd/wpp/DataQuery/

typesoutput = {'lowinc'; 'lmid'; 'umid'};

% Bring in data on life expectancy at birth (for DALYs calculation)
life_expectancy = readtable('life_expectancy.csv','ReadVariableNames',true);

% age groups: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
al = 6;

% Bring in R0, m1, m2, rC, and rep
%load('model_samples.mat')

% New vaccination parameters (same 2000 draws for all countries.)
load('./vax_samp2000dec.mat')

output2 = struct([]);

% what happens in matlab at values of R0 close to 1?
warning('off', 'MATLAB:warn_r14_stucture_assignment')

nruns = 200; %original for Bilcke paper is 2000, set to 30 for test runs

for c=1:ncountry %Loop through countries in test run

country=vimc_countries{c};
%country_name=vimc_country_names{c};

% set z to be the index of whichever country is being modeled
%z = str2num(getenv('OMPI_COMM_WORLD_RANK'))+1;
%z = 48; %index of nigeria in gavi73
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));

%j will map the gavi73 to the lmic 145
%j=gavi73.cn(z); 
j=z; %all_countries.countryn(z); 

tic

% age distribution for this particular country
if agedist_country(j)==1
    dist = low_inc; % change if we change the assumption
elseif agedist_country(j)==2
    dist = lmid_inc; % change if we change the assumption
else
    dist = umid_inc; % change if we change the assumption
end

% Age groups in ODE: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
% Ages for population data: 0-4, 5-9, 10-14, 15-20, 20-25. 
agedist = [dist(1)/5; dist(1)/5; dist(1)*3/5; sum(dist(2:3)); sum(dist(4:5)); 1-sum(dist)];
population = 1e5*agedist;

agepar.mub = [-log(1-mub(agedist_country(j))/1000)/52; zeros(5, 1)]; 
agepar.mu = ([agepar.mub(1); agepar.u(1:(end-1), 1)].*[sum(population); population(1:(end-1),1)])./population(:,1); 
agepar.mu = agepar.mu - agepar.u;

%get population in matrix form for multiplying
vimc_pop_mat(:,:,c) = table2array(input_pop_data(strcmp(input_pop_data.country_code,country),5:105));
%test_pop_mat(:,:,c) = table2array(test5_pop(strcmp(test5_pop.country_code,country),4:104));

%% Calculate multipliers for deaths and Dalys

% same method as Joke uses in R code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antimicrobial Resistance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%proportion of cases AMR (uniform from 0 to 1)
propAMR = rand(nruns,ncountry);

%relative risk of deaths/dalys (uniform from 1 to 3)
rAMR = 1 + 2*rand(nruns,ncountry);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multipliers for Deaths:%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% hospitalized cases = cases*prob hosp
% prob hosp 0·038 (0·038) [0·004-0·249] Common estimate
logit_mean = log(.038/(1-.038));
logit_se = mean([abs(log(.004/(1-.004))-log(.038/(1-.038))),log(.249/(1-.249))-log(.038/(1-.038))])/1.96;
% generate inverse logit estimates and then re-logit to get real values
phosp = normrnd(logit_mean,logit_se,nruns,ncountry);
if strcmp(country,'BDG')
    phosp(:,c) = normrnd(-3.27,1.62,nruns,1);
elseif strcmp(country,'IND')
    phosp(:,c) = normrnd(-2.78,0.37,nruns,1);
elseif strcmp(country,'KEN')
    phosp(:,c) = normrnd(-2.77,1.09,nruns,1);
elseif strcmp(country,'PAK')
    phosp(:,c) = normrnd(-4.21,0.58,nruns,1);
end
%%% Probability of hospitalization %%%%
phosp = exp(phosp)./(1+exp(phosp));

% deaths among hospitalized = hospcases * prob death given hospitalization
pdeathhosp = normrnd(-3.067,0.888,nruns,ncountry);
if country=='BDG'
    pdeathhosp(:,c) = normrnd(-3.45,0.72*2,nruns,1);
elseif country=='ETH'
    pdeathhosp(:,c) = normrnd(-1.99,0.44*2,nruns,1);
elseif country=='CIV'
    pdeathhosp(:,c) = normrnd(-3.24,0.36*2,nruns,1);
elseif country=='LAO'
    pdeathhosp(:,c) = normrnd(-4.27,0.58*2,nruns,1);
elseif country=='ZWE'
    pdeathhosp(:,c) = normrnd(-2.59,0.25*2,nruns,1);
elseif country=='IND'
    pdeathhosp(:,c) = normrnd(-2.80,0.38,nruns,1);
elseif country=='NGA'
    pdeathhosp(:,c) = normrnd(-1.69,0.28,nruns,1);
elseif country=='SEN'
    pdeathhosp(:,c) = normrnd(-3.41,0.51,nruns,1);
end
%%% Probability of Death given hospitalization %%%%
pdeathhosp = exp(pdeathhosp)./(1+exp(pdeathhosp));

% all deaths = hospdeaths/proportion of deaths in hospital
% prop deaths outside hospital is rand uniform between 0 and 75%, so 1-that
% is prop hospital

%%%% Proportion of deaths occuring in hospital %%%%
propdeathhosp = 1-.75*rand(nruns,ncountry);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Multipliers for Dalys:%
%%%%%%%%%%%%%%%%%%%%%%%%%

% years of life affected by disability
% hosp cases * duration of infection hosp *disability weight
% outpatient cases * duration of infection outpatient *disability weight
% cases not seeking med care * doi no med care * disability weight

% prob seeking healthcare = .57, sd .09
% proportion not seeking care = 1- prop seeking care
% prop outpatient = prop seeking care - phosp

% proportion of cases seeking med care 
pmedcare = normrnd(.57,.09,nruns,ncountry);
poutp = pmedcare-phosp;
% find any zeros and replace them with average 
pmedcare(poutp<0) =.57;
phosp(poutp<0) = .038;
poutp(poutp<0) = .57-.038; 
pnomed = 1-pmedcare;

% doi hosp and outpatient = 16 days +-2
doi_hosp = normrnd(16/365.25,2/365.25,nruns,ncountry);
doi_outp = normrnd(16/365.25,2/365.25,nruns,ncountry);
% doi not seeking med care = avg 8 days, 1-100% range of 16 days
doi_nomed = (16/365.25)*rand(nruns,ncountry);
% disability weights: 
% hosp: 0·21 + 0·.04
dwt_hosp = normrnd(.21,.04,nruns,ncountry);
% following bilcke at all, half of runs are simulations with outpatient 
% disability weight for infectious disease, acute, moderate = 0·053 +- 0·012 
% and half for severe = 0·21 + 0·.04
dwt_outp = [normrnd(.053,.012,nruns/2,ncountry);normrnd(.21,.04,nruns/2,ncountry)];
% for those not seeking medical care, half assumed disability weight for
% infectious disease, acute, mild = 0·005 + 0·002
% and half for moderate 0·053 +- 0·012 
dwt_nomed = [normrnd(.005,.002,nruns/2,ncountry);normrnd(.053,.012,nruns/2,ncountry)];

%To get deaths, multiply cases by phosp*pdeathhosp/propdeathosp:
deaths_mult = phosp.*pdeathhosp./propdeathhosp;

%To get ylds, multiply cases by 
% phosp*doihosp*dwthosp + poutp*doioutp*dwtoutp
% +pnomed*doinomed*dwtnomed
yld_mult = phosp.*doi_hosp.*dwt_hosp+poutp.*doi_outp.*dwt_outp+pnomed.*doi_nomed.*dwt_nomed;

% years of life lost = deaths*life expectancy - age at death
%set up vector of years of life lost for each age
yr_of_birth=zeros(101,101);
for y=1:101
    for a=1:101
        yr_of_birth(y,a)=2000+y-a;
    end
end

leb_mat=zeros(101,101);
leb_country=life_expectancy(strcmp(life_expectancy.country_code,country),:);
for y=1:101
    for a=1:101
        if yr_of_birth(y,a)<=1950 
            leb_mat(y,a)=leb_country.value(1);
        elseif yr_of_birth(y,a)>2099
            leb_mat(y,a)=leb_country.value(end);
        else
            leb_mat(y,a)=leb_country.value(find(leb_country.year==yr_of_birth(y,a)));
        end
    end
end

%yll_vec = repmat([leb-[0:69],zeros(1,31)],101,1);
yll_vec = zeros(101,101);
for a=1:101
    yll_vec(:,a) = max(0,leb_mat(:,a)-a+1);
end

%% Generate stochastic simulations

tic  

for i=1:nruns 
    
estpar.mult1 = double(m1_sample(i,1)); 
estpar.mult2 = double(m2_sample(i,1)); 
estpar.r = 1; 
estpar.rC = double(rC_sample(i,1)); 
coeff = [estpar.mult1,estpar.mult1,estpar.mult2,1,1,1]; 
estpar.R0 = double(R0samples(i,z)); 
estpar.beta = estpar.R0*coeff'.*(agepar.mub(1) + params.delta)./((1+params.delta*(agedist'*agepar.theta)*estpar.rC./agepar.mub(1))*(coeff*agedist));
params.epsilon = 0;
h2o=unifrnd(-0.025,0);
    
vacpar.omega_v = -log(1-vax_samp.omega(i))/52; % convert to WEEKS
vacpar.veff = vax_samp.nu(i);

    % Simulate equilibrium
    % Run to equilibrium only once and then evaluate both unc and con and all interventions with this.
    options = odeset('RelTol',1e-6);
    
    % How long to reach equilibrium
    if estpar.beta(4)<0.2
    tspan=7000*52+1;
    elseif estpar.beta(4)<0.3
    tspan=3000*52+1;
    else
    tspan=200*52+1;
    end
    
    pop0 = zeros(al*18,1);
    pop0(1:al) = population; % putting everyone in susceptible 2???
    % throw an infected in one of the four categories to start the simulation
    ageptzero = 3; %randsample(1:al, 1, true, pop0(1:al)/sum(pop0(1:al)));
    pop0(ageptzero) = pop0(ageptzero)-100;
    pop0(al+ageptzero) = 100;
    
    rvnone = zeros(tspan, 6);
    massvacnone = zeros(tspan, 6);
    
    %estpar.beta_h2o = repmat(estpar.beta, 1, 1+52*10).*repmat(1:-0.25/length(1:(52*10)):0.75, al, 1);
    %estpar.beta_h2o = repmat(estpar.beta, 1, 1+52*int_length).*repmat(1:-0.25/length(1:(52*int_length)):0.75, al, 1);
    estpar.beta_h2o = repmat(estpar.beta, 1, 52*101+1).*repmat(exp(h2o*(-15.5:(1/52):85.5)), al, 1);
    
    %Burn-in period to reach equilibrium
    [t, pop] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata(t, pop, estpar.beta_h2o(:,1), agepar.mub, agepar.u, ...
        agepar.mu, params.omega, vacpar.omega_v, estpar.r, estpar.rC, params.alpha, params.delta, ...
        agepar.theta, agepar.theta2, params.epsilon, rvnone, vacpar.veff, massvacnone, al, population), ...
        1:tspan, pop0, options);
    pop0eqm = pop(end, :);
    %save initial population
    pre_pop = pop;
    totpop = sum(pop,2); %get total population
    %get age-specific incidence
    pre_inctime_unvacc = diff(pop(:,(16*al+1):(17*al)));
     
    
    % For VIMC, not doing different coverage scenarios, just given coverage
    cov_est = vimc_cov;
    
    % coverage - this time it is going to be from the starting point until
    % 30 years after the starting point...
    tmp = repmat(cov_est.coverage(strcmp(cov_est.country_code, iso_lmic.countryiso(j)) & strcmp(cov_est.activity_type, 'routine')), 1, 52);
    tmp(isnan(tmp(:,1)), :) = [];
    
    pre_int_length = find(tmp,1); % length of time before vaccine roll-out
    int_length = 101-pre_int_length; %height(rfpNGAcov); %number of years post-intervention
    
    tmp = tmp';
    tmp = tmp(:);

    if length(tmp)>1
        tmp1 = [0; tmp];
    else
        tmp1 = zeros(5253,1);
    end
    
    v1 = [zeros(length(tmp1),1), tmp1, zeros(length(tmp1),4)];
    % only the second age group gets vaccination, hence 0's for every other
    % age group
    
    tmp2 = cov_est.coverage(strcmp(cov_est.country_code, iso_lmic.countryiso(j)) & strcmp(cov_est.activity_type, 'campaign'));    
    
    tmp2(isnan(tmp2)) = [];
    camp_ages = zeros(length(tmp1), 1);
    
    % this is coded to tolerate the scenario when campaigns are done over
    % the span of several years
    for k=1:length(tmp2)
        camp_ages((52*(pre_int_length+k-1)+1):(52*(pre_int_length+k-1)+4)) = min((1-(1-tmp2(k)).^.25),.62);
    end

    % age groups: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
    massvac{1,1} = zeros(length(camp_ages), 6); % No vaccination
    massvac{2,1} = [zeros(length(camp_ages), 1), repmat(camp_ages, 1, 3), zeros(length(camp_ages), 2)]; % Campaign only
    massvac{3,1} = [zeros(length(camp_ages), 1), repmat(camp_ages, 1, 3), zeros(length(camp_ages), 2)]; % Routine + campaign
    
    for int=1:3 % cycle through interventions
        if int<3
            rv = zeros(size(v1,1), size(v1,2));
        else
            rv=v1;
        end
        
        if int==1 || (int>1 && length(tmp)>1)  % Don't need to simulate the vaccination scenarios for countries with coverage=0 
        
    [t, pop] = ode45(@(t, pop) fn_SIR_vaccine_vnv_strata_h2o(t, pop, estpar.beta_h2o, agepar.mub, agepar.u, ...
        agepar.mu, params.omega, vacpar.omega_v, estpar.r, estpar.rC, params.alpha, params.delta, ...
        agepar.theta, agepar.theta2, params.epsilon, rv, vacpar.veff, massvac{int, 1}, al, population), ...
        1:(length(tmp1)), pop0eqm, options);
    
        end
    
    % chronic = pop(5:end, (5*al+1):(6*al)) + pop(:, (13*al+1):(14*al));    
    inctime_unvacc = diff(pop(:,(16*al+1):(17*al))); %Incidence among unvaccinated
    inctime_vacc = diff(pop(:,(17*al+1):(18*al))); %Incidence among vaccinated
    %get total incidence by age, weekly
    inctime_tot = inctime_unvacc+inctime_vacc;
    %get total incidence by age, year
    inctime_age = squeeze(single(sum(reshape(inctime_tot, 52, pre_int_length+int_length, al),1)));
    %get population size
    tmp_pop=reshape(pop(:,1:al*14),length(pop),al,14); %population by age and disease state
    age_pop = sum(tmp_pop,3); %sum to get total population by age
    age_pop = age_pop(26:52:end,:); %get midpoint population of each year for division
    
    %get age-specific incidence rate per individual
    inc_rate_age = inctime_age./age_pop;
    %get age-specific incidence rate per 100K
    inc_rate_age_100K=inc_rate_age*1e5;
    
    %calculate individual year age-specific incidence
    % distribute out incidence levels to midpoint of age group
    %try endpoint
    %Age groups in ODE: 0-9m, 9m-2y, 2y-4y, 5y-14y, 15y-25y, 25y+
    x_age = [0 .375 1.375 3.5 10 20 62.5 100];
    
    % interpolate individual year age-specific incidence rate
    interp_inc=pchip(x_age,[inc_rate_age(:,1) inc_rate_age inc_rate_age(:,6)],0:100);
    

    rvac = diff(sum(pop(:,(14*al+1):(15*al)),2));
    cvac = diff(sum(pop(:,(15*al+1):(16*al)),2));
    
    % Store: 
    % output2{i,cov}{int,1} = single(1);
    output2{i,c}{int,1}.interp_inc=interp_inc; %interpolated individual year age-specific incidence rate
    output2{i,c}{int,1}.inc_rate_age=inc_rate_age; %age-group age-specific incidence rate
    output2{i,c}{int,1}.inc_rate_age_100K=inc_rate_age_100K; %age-group age-specific incidence rate
    output2{i,c}{int,1}.inc_rate_100K=sum(inctime_age,2)./sum(age_pop,2)*1e5;

    output2{i,c}{int,1}.inctime_age=inctime_age; %age-specific incidence by year, vacc+unvacc
    output2{i,c}{int,1}.age_pop=age_pop;  %age-specific population size over time
%     for age=1:6
%         %get age-specific incidence
%         output2{i,cov}{int,1}.inctime_unvacc(:,age) = single(sum(reshape(inctime_unvacc(:,age), 52, 10), 1));
%         output2{i,cov}{int,1}.inctime_vacc(:,age) = single(sum(reshape(inctime_vacc(:,age), 52, 10),1));
%     
%     end
   
    output2{i,c}{int,1}.rvac = single(sum(reshape(rvac, 52, 101), 1));
    output2{i,c}{int,1}.cvac = single(sum(reshape(cvac, 52, 101), 1));
    output2{i,c}{int,1}.mult1 = single(estpar.mult1); 
    output2{i,c}{int,1}.mult2 = single(estpar.mult2);
    output2{i,c}{int,1}.rC = single(estpar.rC);
    output2{i,c}{int,1}.beta = single(estpar.beta);
    output2{i,c}{int,1}.R0 = single(estpar.R0);
    output2{i,c}{int,1}.omega_v = single(vacpar.omega_v);
    output2{i,c}{int,1}.veff = single(vacpar.veff);
    output2{i,c}{int,1}.h2o = h2o;
   
    %% Calculate cases, deaths and dalys 
    output2{i,c}{int,1}.cases = round(repsamples(i,z)*output2{i,c}{int,1}.interp_inc.*vimc_pop_mat(:,:,c)');
    output2{i,c}{int,1}.deaths = round((1-propAMR(i)+rAMR(i)*propAMR(i))*deaths_mult(i)*output2{i,c}{int,1}.cases);
    output2{i,c}{int,1}.dalys = (1-propAMR(i)+rAMR(i)*propAMR(i))*(yld_mult(i)*output2{i,c}{int,1}.cases+...
        deaths_mult(i)*output2{i,c}{int,1}.cases.*yll_vec);
    end
end

end

toc

%filename = ['./output2_testcountry', sprintf('%02.0f', j),'_h2o.mat'];
save('output2_vimc_full.mat'); 
%end


%% Generate figure of interpolated incidence for 2018

test_countries = {'COD','ETH','IND','NGA','PAK'};

set(0,'defaultAxesFontSize',10)

figure('units','normalized','outerposition',[0 0 .5 .5])

interp_mat = zeros(101,nruns,5);
age_inc_mat = zeros(6,nruns,5);
x_age_plot = [.375 1.375 3.5 10 20 62.5];

for c=1:5
    for i = 1:nruns
        interp_mat(:,i,c) = output2{i,c}{1,1}.interp_inc(19,:);
        age_inc_mat(:,i,c) = output2{i,c}{1,1}.inc_rate_age(19,:);
    end
%generate means and error bars for original age-specific incidence and
%interpolated incidenc
mean_interp = 1e5*mean(interp_mat(:,:,c),2);
neg_err_interp = 1e5*(mean(interp_mat(:,:,c),2)-prctile(interp_mat(:,:,c),2.5,2));
pos_err_interp = -1e5*(mean(interp_mat(:,:,c),2)+prctile(interp_mat(:,:,c),97.5,2));
mean_age_inc = 1e5*mean(age_inc_mat(:,:,c),2);
neg_err_age_inc = 1e5*(mean(age_inc_mat(:,:,c),2)-prctile(age_inc_mat(:,:,c),2.5,2));
pos_err_age_inc = -1e5*(mean(age_inc_mat(:,:,c),2)+prctile(age_inc_mat(:,:,c),97.5,2));

subplot(2,3,c)
hold on
errorbar([0:100],mean_interp,neg_err_interp,pos_err_interp,'Color',[.7 .7 .7])
errorbar(x_age_plot,mean_age_inc',neg_err_age_inc,pos_err_age_inc,'ro','MarkerSize',5,'MarkerFaceColor','r');
plot([0:100],mean_interp,'b','LineWidth',2);
xlim([-1 101])
ylabel('Incidence per 100K py')
title({'Incidence Rate by Age'; test_countries{c}})
%title({'Incidence Rate by Age'; vimc_countries{c}})
end
legend('Interpolated','Model Original')


%% Plot age-specific and overall incidence per 100K

age_grps = {'0-9m', '9m-2y', '2y-4y', '5y-14y', '15y-25y', '25y+'};

figure

for c=1:5
country=test_countries{c};
%country=vimc_countries{c};
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));
age_inc_baseline = [];
age_inc_RC15= [];
age_inc_camp= [];

for i=1:size(output2,1)
    age_inc_baseline=cat(3,age_inc_baseline, repsamples(i,z)*output2{i,c}{1,1}.inc_rate_age_100K);
    age_inc_camp=cat(3,age_inc_camp, repsamples(i,z)*output2{i,c}{2,1}.inc_rate_age_100K);         
    age_inc_RC15=cat(3,age_inc_RC15, repsamples(i,z)*output2{i,c}{3,1}.inc_rate_age_100K);
end

for j = 1:6
    subplot(6,5,5*(j-1)+c)
    hold on
    plot(mean(age_inc_baseline(:,j,:),3),'b','LineWidth',1);
    plot(median(age_inc_baseline(:,j,:),3),'r','LineWidth',1);
    plot(mean(age_inc_RC15(:,j,:),3),'--b','LineWidth',1);
    plot(median(age_inc_RC15(:,j,:),3),'--r','LineWidth',1);
    plot(mean(age_inc_camp(:,j,:),3),':b','LineWidth',1);
    plot(median(age_inc_camp(:,j,:),3),':r','LineWidth',1);
    for i=1:size(output2,1)
        plot(repsamples(i,z)*output2{i,c}{1,1}.inc_rate_age_100K(:,j),'Color',[.7 .7 .7])
    end
    plot(mean(age_inc_baseline(:,j,:),3),'b','LineWidth',1);
    plot(median(age_inc_baseline(:,j,:),3),'r','LineWidth',1);
    plot(mean(age_inc_RC15(:,j,:),3),'--b','LineWidth',1);
    plot(median(age_inc_RC15(:,j,:),3),'--r','LineWidth',1);
    plot(mean(age_inc_camp(:,j,:),3),':b','LineWidth',1);
    plot(median(age_inc_camp(:,j,:),3),':r','LineWidth',1);
    set(gca,'XLim',[0 100],'XTick',0:50:100,'XTickLabel',2000:50:2100)
    if j==3
    ylabel('Incidence per 100K PY')
    end
    %title([age_grps{j}, vimc_countries{c}])
    title([age_grps{j}, test_countries{c}])
end
end
legend('Baseline Mean','Baseline Median','RoutineCamp Mean','RoutineCamp Median','CampOnly Mean','CampOnly Median')
    
%%

inc_baseline = zeros(101,size(output2,1));
inc_camp= zeros(101,size(output2,1));
inc_RC15= zeros(101,size(output2,1));

figure
for c=1:5
country=test_countries{c};
%country=vimc_countries{c};
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));

for i=1:size(output2,1)
    inc_baseline(:,i) = repsamples(i,z)*output2{i,c}{1,1}.inc_rate_100K;
    inc_camp(:,i) = repsamples(i,z)*output2{i,c}{2,1}.inc_rate_100K;
    inc_RC15(:,i) = repsamples(i,z)*output2{i,c}{3,1}.inc_rate_100K;      
end

subplot(3,2,c)
hold on
plot(mean(inc_baseline,2),'b','LineWidth',2);
plot(median(inc_baseline,2),'r','LineWidth',2);
plot(mean(inc_RC15,2),'--b','LineWidth',2);
plot(median(inc_RC15,2),'--r','LineWidth',2);
plot(mean(inc_camp,2),':b','LineWidth',2);
plot(median(inc_camp,2),':r','LineWidth',2);
for i=1:size(output2,1)
    plot(repsamples(i,z)*output2{i,c}{1,1}.inc_rate_100K,'Color',[.7 .7 .7])
end
plot(mean(inc_baseline,2),'b','LineWidth',2);
plot(median(inc_baseline,2),'r','LineWidth',2);
plot(mean(inc_RC15,2),'--b','LineWidth',2);
plot(median(inc_RC15,2),'--r','LineWidth',2);
plot(mean(inc_camp,2),':b','LineWidth',2);
plot(median(inc_camp,2),':r','LineWidth',2);
if c==1
legend('Baseline Mean','Baseline Median','RoutineCamp Mean','RoutineCamp Median','CampOnly Mean','CampOnly Median')
end
set(gca,'XLim',[0 100],'XTick',0:20:100,'XTickLabel',2000:20:21000)
title({'Incidence per 100K PY, Total'; vimc_countries{c}})
end

%% Check to make sure overall incidence and aggregated overall incidence 
% from age groups match roughly, baseline only

%get total population over time

figure
x=[16 23 36 57 60]; % VIMC index of test countries

for c=1:5
country=test_countries{c};
%country=vimc_countries{c};
test_tot_pop = sum(vimc_pop_mat(:,:,x(c)),1)';
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));

tot_inc = zeros(101,nruns);
age_ag_inc = zeros(101,nruns);
for i = 1:nruns
    tot_inc(:,i) = (repsamples(i,z)*output2{i,x(c)}{1,1}.inc_rate_100K/1e5).*test_tot_pop;
    age_ag_inc(:,i)= repsamples(i,z)*sum(output2{i,x(c)}{1,1}.interp_inc'.*vimc_pop_mat(:,:,x(c)),1)';
end
subplot(2,3,c)
hold on
plot(mean(tot_inc,2),'b')
plot(mean(age_ag_inc,2),'r')
plot(prctile(tot_inc,2.5,2),'--b')
plot(prctile(tot_inc,97.5,2),'--b')
plot(prctile(age_ag_inc,2.5,2),'--r')
plot(prctile(age_ag_inc,97.5,2),'--r')
%title(vimc_countries{c})
set(gca,'XLim',[0 102],'XTick',1:20:101,'XTickLabel',2000:20:2100,'XTickLabelRotation',45)
ylabel('Incidence per 100K py')
title(test_countries{c})
end
legend('Total Incidence','Aggregated Age Incidence')

%% Plot the country-specific age distribution in 2012 vs the modeled age distribution

figure
x=[16 23 36 57 60]; % VIMC index of test countries

for c=1:5
country=test_countries{c};
%country=vimc_countries{c};
test_pop_agedist = [sum(vimc_pop_mat(1:5,13,x(c)),1) sum(vimc_pop_mat(6:10,13,x(c)),1) sum(vimc_pop_mat(11:15,13,x(c)),1) sum(vimc_pop_mat(15:20,13,x(c)),1) sum(vimc_pop_mat(20:24,13,x(c)),1) sum(vimc_pop_mat(25:end,13,x(c)),1)]/sum(vimc_pop_mat(:,13,x(c)),1);
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));

if agedist_country(z)==1
    dist = [low_inc 1-sum(low_inc)]; % change if we change the assumption
elseif agedist_country(z)==2
    dist = [lmid_inc 1-sum(lmid_inc)]; % change if we change the assumption
else
    dist = [umid_inc 1-sum(umid_inc)]; % change if we change the assumption
end

subplot(2,3,c)
bar([dist' test_pop_agedist'])
set(gca,'XTickLabel',{'0-4y','5-9y','10-14y','15-20y','20-24y','25+y'},'XTickLabelRotation',45)
xlabel('Age group')
ylabel('Prop. of population')
title(test_countries{c})
end
legend('Modeled population','Country-specific age dist (2012)')

%%
tot_cases_novacc=zeros(101,nruns);
tot_cases_RC15=zeros(101,nruns);
tot_cases_camp=zeros(101,nruns);
tot_cases_all_novacc=zeros(101,nruns);
tot_cases_all_RC15=zeros(101,nruns);
tot_cases_all_camp=zeros(101,nruns);

figure
for c=1:5
country=test_countries{c};
%country=vimc_countries{c};
%z = find(strcmp(gavi73.ISO3,country));
z = find(strcmp(all_countries.ISO3,country));

for i=1:nruns
    tot_cases_novacc(:,i)=sum(output2{i,c}{1,1}.cases,2);
    tot_cases_camp(:,i)=sum(output2{i,c}{2,1}.cases,2);
    tot_cases_RC15(:,i)=sum(output2{i,c}{3,1}.cases,2);
end
tot_cases_all_novacc=tot_cases_all_novacc+tot_cases_novacc;
tot_cases_all_camp=tot_cases_all_camp+tot_cases_camp;
tot_cases_all_RC15=tot_cases_all_RC15+tot_cases_RC15;

subplot(3,2,c)
hold on
plot(mean(tot_cases_novacc,2),'b','LineWidth',2);
plot(mean(tot_cases_RC15,2),'r','LineWidth',2);
plot(mean(tot_cases_camp,2),'g','LineWidth',2);
for i=1:size(output2,1)
    plot(sum(output2{i,c}{1,1}.cases,2),'Color',[.7 .7 .7])
end
plot(mean(tot_cases_novacc,2),'b','LineWidth',2);
plot(mean(tot_cases_RC15,2),'r','LineWidth',2);
plot(mean(tot_cases_camp,2),'g','LineWidth',2);
if c==1
legend('NoVacc Mean','RoutineCamp Mean','CampOnly Mean')
end
set(gca,'XLim',[0 100],'XTick',0:20:100,'XTickLabel',2000:20:21000)
%title({'TF Cases per year, Total'; vimc_countries{c}})
title({'TF Cases per year, Total'; test_countries{c}})
end

subplot(3,2,6)
hold on
plot(mean(tot_cases_all_novacc,2),'b','LineWidth',2);
plot(mean(tot_cases_all_RC15,2),'r','LineWidth',2);
plot(mean(tot_cases_all_camp,2),'g','LineWidth',2);
if c==1
legend('NoVacc Mean','RoutineCamp Mean','CampOnly Mean')
end
set(gca,'XLim',[0 100],'XTick',0:20:100,'XTickLabel',2000:20:21000)
title({'TF Cases per year, Total'; 'All countries'})

