%% Generate output tables formatted for VIMC

% For stochastic outputs, one for baseline (no vacc), one for camp 15 (default):
%for each run id, need year in order, then age group, cohort size, cases, dalys
%stochastic_output_novacc = table; 

% For central output: same as stochastic, but mean values
central_output_novacc = table; 

% read in disease, years, age, country, cohort size, and then set up default tables as copies
% stochastic output tables
stochastic_output_novacc.disease(1:nruns*ncountry*101*101,1) = {'Typhoid'};

stochastic_output_campaign = stochastic_output_novacc;
stochastic_output_camproutine = stochastic_output_novacc;
central_output_campaign = central_output_novacc;
central_output_camproutine = central_output_novacc;

for c=1:ncountry
    
country=vimc_countries{c}; %test_countries{c};
country_name=vimc_country_names{c};   
    
stochastic_output_novacc.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_novacc.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_novacc.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_novacc.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_novacc.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_novacc.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_campaign.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_campaign.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_campaign.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_campaign.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_campaign.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_campaign.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

stochastic_output_camproutine.run_id = repmat(reshape(repmat(1:nruns,101*101,1),nruns*101*101,1),ncountry,1);
stochastic_output_camproutine.year =reshape(repmat(2000:2100,1,nruns*ncountry*101),1,nruns*ncountry*101*101)';
stochastic_output_camproutine.age = reshape(repmat(0:100,101,nruns*ncountry),nruns*ncountry*101*101,1);
stochastic_output_camproutine.country(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country};
stochastic_output_camproutine.country_name(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = {country_name};
stochastic_output_camproutine.cohort_size(nruns*101*101*(c-1)+1:nruns*101*101*c,1) = repmat(reshape(vimc_pop_mat(:,:,c)',101*101,1),nruns,1);

%central output tables
central_output_novacc.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_novacc.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_novacc.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_novacc.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_novacc.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_novacc.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_campaign.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_campaign.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_campaign.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_campaign.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_campaign.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_campaign.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);

central_output_camproutine.disease(1:ncountry*101*101,1) = {'Typhoid'};
central_output_camproutine.year = reshape(repmat(2000:2100,1,ncountry*101),ncountry*101*101,1);
central_output_camproutine.age = reshape(repmat(0:100,101,ncountry),ncountry*101*101,1);
central_output_camproutine.country(101*101*(c-1)+1:101*101*c,1) = {country};
central_output_camproutine.country_name(101*101*(c-1)+1:101*101*c,1) = {country_name};
central_output_camproutine.cohort_size(101*101*(c-1)+1:101*101*c,1) = reshape(vimc_pop_mat(:,:,c)',101*101,1);


%for c=1:ncountry
for i =1:size(output2,1)  
    % Read in cases, deaths, dalys for run_id ==1
    % no vacc
    %stochastic_output_novacc.cases(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.cases,101*101,1);
    %stochastic_output_novacc.deaths(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.deaths,101*101,1);
    %stochastic_output_novacc.dalys(logical((stochastic_output_novacc.run_id==i).*strcmp(stochastic_output_novacc.country,country)),1) = reshape(output2{i,c}{1,1}.dalys,101*101,1);
    stochastic_output_novacc.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{1,1}.cases,101*101,1);
    stochastic_output_novacc.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{1,1}.deaths,101*101,1);
    stochastic_output_novacc.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{1,1}.dalys,101*101,1);

    % Campaign only
    stochastic_output_campaign.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{2,1}.cases,101*101,1);
    stochastic_output_campaign.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{2,1}.deaths,101*101,1);
    stochastic_output_campaign.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{2,1}.dalys,101*101,1);

    % Campaign & Routine
    stochastic_output_camproutine.cases(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{3,1}.cases,101*101,1);
    stochastic_output_camproutine.deaths(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{3,1}.deaths,101*101,1);
    stochastic_output_camproutine.dalys(nruns*101*101*(c-1)+101*101*(i-1)+1:nruns*101*101*(c-1)+101*101*i,1) = reshape(output2{i,c}{3,1}.dalys,101*101,1);

end

% generate means for central output
central_output_novacc.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_novacc.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_novacc.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_novacc.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_novacc.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_novacc.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

central_output_campaign.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_campaign.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_campaign.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_campaign.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

central_output_camproutine.cases(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine.cases(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_camproutine.deaths(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine.deaths(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);
central_output_camproutine.dalys(101*101*(c-1)+1:101*101*c,1) = reshape(mean(reshape(stochastic_output_camproutine.dalys(nruns*101*101*(c-1)+1:nruns*101*101*c,1),101,101,nruns),3),101*101,1);

end

%mean_total_cases_novacc = median(squeeze(sum(reshape(stochastic_output_novacc.cases,101,101,nruns),2)),2);
%mean_total_deaths_novacc = median(squeeze(sum(reshape(stochastic_output_novacc.deaths,101,101,nruns),2)),2);
%mean_total_dalys_novacc = median(squeeze(sum(reshape(stochastic_output_novacc.dalys,101,101,nruns),2)),2);

%mean_total_cases_default = median(squeeze(sum(reshape(stochastic_output_default.cases,101,101,nruns),2)),2);
%mean_total_deaths_default = median(squeeze(sum(reshape(stochastic_output_default.deaths,101,101,nruns),2)),2);
%mean_total_dalys_default = median(squeeze(sum(reshape(stochastic_output_default.dalys,101,101,nruns),2)),2);

%total_output_novacc = table([2000:2100]',mean_total_cases_novacc, mean_total_deaths_novacc, mean_total_dalys_novacc,...
%    'VariableName',{'Year','MeanTotalCases','MeanTotalDeaths','MeanTotalDalys'});
%total_output_default = table([2000:2100]',mean_total_cases_default, mean_total_deaths_default, mean_total_dalys_default,...
%    'VariableName',{'Year','MeanTotalCases','MeanTotalDeaths','MeanTotalDalys'});

writetable(stochastic_output_novacc,'stochastic_output_TF-Yale-Pitzer_novacc.csv')
writetable(stochastic_output_campaign,'stochastic_output_TF-Yale-Pitzer_campaign.csv')
writetable(stochastic_output_camproutine,'stochastic_output_TF-Yale-Pitzer_camproutine.csv')

writetable(central_output_novacc,'central_output_TF-Yale-Pitzer_novacc.csv')
writetable(central_output_campaign,'central_output_TF-Yale-Pitzer_campaign.csv')
writetable(central_output_camproutine,'central_output_TF-Yale-Pitzer_camproutine.csv')

%writetable(total_output_novacc,'total_output_novacc.csv')
%writetable(total_output_default,'total_output_default.csv')

%% Compile parameter file to export

param_table = table; %'Size',[nruns,27],...

param_table.run_id = (1:nruns)';
for i = 1:nruns
    for c=1:ncountry
        z = find(strcmp(all_countries.ISO3,vimc_countries{c}));
        eval(strcat('param_table.R0_',vimc_countries{c},'(i)=','R0samples(i,z);'));
        eval(strcat('param_table.rep_',vimc_countries{c},'(i)=','repsamples(i,z);'));
    end
    %read in parameters that are in output2 structure
    %param_table.R0(i) = output2{i,1}{1,1}.R0;
    param_table.mult1(i) = output2{i,c}{1,1}.mult1;
    param_table.mult2(i) = output2{i,c}{1,1}.mult2;
    param_table.rC(i) = output2{i,c}{1,1}.rC;
    param_table.h2o(i) = output2{i,c}{1,1}.h2o;
    param_table.veff(i) = output2{i,c}{1,1}.veff;
    param_table.omega_v(i) = vax_samp.omega(i);
    %param_table.delta(i) = params.delta;
    %param_table.alpha(i) = params.alpha;
    %param_table.omega(i) = params.omega;
    %param_table.mub(i) = mub(agedist_country(gavi73.cn(z))); %birth rate
    %param_table.mub(i) = mub(agedist_country(all_countries.countryn(z))); %birth rate
    %param_table.theta_y(i) = agepar.theta(1);
    %param_table.theta_o(i) = agepar.theta(end);
    %param_table.r(i) = 1;
    %param_table.leb(i)= leb;
    for str2 = {'propAMR','rAMR','phosp','pmedcare','pdeathhosp','propdeathhosp','doi_hosp','doi_outp',...
    'doi_nomed','dwt_hosp','dwt_outp','dwt_nomed'}
        eval(strcat('param_table.',str2{:},'(i)=',str2{:},'(i);'));
    end
    
end
writetable(param_table,'TF-Yale-Pitzer_stochastic_parameters.csv');