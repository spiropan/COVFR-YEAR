% COVID-19_vaccinations_sorted downloaded from
% https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc

% Deaths by sex and age downloaded from
% https://data.cdc.gov/NCHS/Provisional-COVID-19-Deaths-by-Sex-and-Age/9bhg-hcku

% 2021 total state populations obtained from:
% https://fred.stlouisfed.org/release/tables?rid=118&eid=259194&od=2021-01-01#
% 2021 State populations over 65 obtained from
% https://www.consumeraffairs.com/homeowners/elderly-population-by-state.html#ratio-by-state

%% RUN THE YEARLY REGRESSIONS
clear
kill_figures
death_types={'COVID','total'};
out.death_types=death_types;
% If below is set to 1, will use robust regression when fitting line
robust=1; 
% If below is set to 0, the code will use y21 vax to predict change in deaths in y22,
% If it is set to 1, the code will use y22 vax to predict change in deaths in y23
do_y23=0; 
% If below is set to 0, the code will do analysis for all ages, if set to
% 1, it will run it only for ages 65 and over
do_over65=1;
% If this is above zero, then will bootstrap to compute 95% CI
iters=2000; 
% If below is set to 1 will regress (dth2-dth1)/dth1 on vax./pop 
perc_reg=0; 
% Percentile outliers to remove, if any. If 0, no outliers will be removed.
% If i.e. 10, then residuals below 10% and above 90% are removed. 
outper=10;

for dt=1:length(death_types)
    disp(['Working on death type: ' death_types{dt}])
    fname=['Table-YEARS-' death_types{dt} '.csv'];
    tname='Tables';
    tab=readtable(['./' tname '/' fname]);
    
    % extract vaccine doses in 2021 and 2022
    if do_over65
        vax_21=tab.(['Y21_Administered_65Plus']);
        vax_22=tab.(['Y22_Administered_65Plus']);
    else
        vax_21=tab.(['Y21_Administered']);
        vax_22=tab.(['Y22_Administered']); 
    end
    
    % extract death counts
    if do_y23
        vax_orig=vax_22;
        if do_over65
            dth1_orig=tab.(['Y22_65_74'])+tab.(['Y22_75_84'])+tab.(['Y22_85_plus']);
            dth2_orig=tab.(['Y23_65_74'])+tab.(['Y23_75_84'])+tab.(['Y23_85_plus']);
        else
            dth1_orig=tab.(['Y22_All']);
            dth2_orig=tab.(['Y23_All']);
        end
        % Scale to account for fact the Y23 goes only through September 27,
        % 2023
        dth2_orig=365/270*dth2_orig;
    else
        vax_orig=vax_21;
        if do_over65
            dth1_orig=tab.(['Y21_65_74'])+tab.(['Y21_75_84'])+tab.(['Y21_85_plus']);
            dth2_orig=tab.(['Y22_65_74'])+tab.(['Y22_75_84'])+tab.(['Y22_85_plus']);
        else
            dth1_orig=tab.(['Y21_All']);
            dth2_orig=tab.(['Y22_All']);
        end
    end
    % extract age-stratified population
    if do_over65
        pop_orig=tab.(['Y21_pop_65Plus']);
    else
        pop_orig=tab.(['Y21_pop_All']);
    end
    
    % Exclude Nans  etc.
    exclude = @(val) val <= 0 | isnan(val);
    rm_idx=exclude(vax_orig) | exclude(dth1_orig) | exclude(dth2_orig) | exclude(pop_orig);
    
    vax_orig(rm_idx)=[]; dth1_orig(rm_idx)=[]; dth2_orig(rm_idx)=[]; pop_orig(rm_idx)=[];    
    %vax_21_orig(rm_idx)=[]; vax_22_orig(rm_idx)=[];
    
    % extract state names
    states_orig=tab.(['State']); states_orig(rm_idx)=[];
    
    % Compute bootstrap estimates
    for i=1:iters+1
        disp(['working on iteration: ' int2str(i)])
        
        % Boostrap if i > 1
        if i==1
            pind=[1:length(dth1_orig)];
        else
            pind = randsample([1:length(dth1_orig)],length(dth1_orig),true);
        end
        
        % Replace for bootstrap indeces if i > 1
        vax=vax_orig(pind); dth1=dth1_orig(pind); dth2=dth2_orig(pind); pop=pop_orig(pind); states=states_orig(pind);
        %vax_21=vax_21_orig(pind); vax_22=vax_22_orig(pind);
        
        % Take the logs
        vax_log=log(vax); dth1_log=log(dth1); dth2_log=log(dth2); pop_log=log(pop);
        
        % Remove outliers to improve VFR estimates
        if perc_reg
            [~,tmp]=robustfit([vax./pop],[(dth2-dth1)./dth1]);
        else
            [~,tmp]=robustfit([dth1_log,vax_log],[dth2_log]);
        end
        [~,I]=rmoutliers(tmp.resid,'percentiles',[outper 100-outper]);
        vax(I)=[]; dth1(I)=[]; dth2(I)=[]; pop(I)=[]; states(I)=[];
        %vax_21(I)=[]; vax_22(I)=[];
        
        % Take the logs
        vax_log=log(vax); dth1_log=log(dth1); dth2_log=log(dth2); pop_log=log(pop);
        
        if perc_reg
            X=vax./pop; Y=(dth2-dth1)./dth1;
        else
            X=[dth1_log vax_log]; Y=dth2_log;
        end
        
        % Estimate the MLR regression models
        if robust  
            [b,stats]=robustfit([X],[Y]);
            [bres,resid]=robustfit([dth1_log],[dth2_log]);
        else
            [b,dev,stats]=glmfit([X],[Y]);
            [~,~,resid]=glmfit([dth1_log],[dth2_log]);
        end
        
        % Store the beta estimates and p-values for the vaccination term
        out.beta(i,dt)=b(end);
        out.pval(i,dt)=stats.p(end);
        
        % compute correlations between residuals and population and
        % percent vaccinated in 2021
        [out.r_pop_resid(i,dt),out.p_pop_resid(i,dt)]=corr(pop,resid.resid,'type','spearman');
        [out.r_per_resid(i,dt),out.p_per_resid(i,dt)]=corr(vax./pop,resid.resid,'type','spearman');
        
        % compute correlations between yearly percent change in deaths and
        % 2021 percent vaccinated
        [out.r_per_per(i,dt),out.p_per_per(i,dt)]=corr(vax./pop,(dth2-dth1)./dth1,'type','spearman');
        
        % Make the plots if i=1
        if i==1
            % Create vectors for regression plot
            if isstruct(resid)
                xvec=vax_log; yvec=resid.resid; %         
            else
                xvec=NaN; yvec=NaN;
            end
            
            % Create vectors for correlation plots
            xvec2=vax./pop; yvec2=(dth2-dth1)./dth1;  
            
            if do_y23
                yr1='2022'; yr2='2023';
            else
                yr1='2021'; yr2='2022';
            end
            
            if do_over65
                ages='65 plus';
            else
                ages='all ages';
            end
            
            % Make figures for the regression
            if ~perc_reg
                figure
                fname1  = ['Fig-regression-' death_types{dt} '-YEAR-' yr2 '-' strrep(ages,' ','-')];
                xlab   = ['Log vaccine doses administered in ' yr1 ' for ' ages];
                ylab   = ['Residual deaths: log(2022 vs 2021), ' ages];
                ftitle = ['Doses vs. ' death_types{dt} ' deaths: tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];
                
                try
                    plot_fit(xvec,yvec,robust,'k.','r',states);
                catch
                    plot_fit(xvec,yvec,0,'k.','r',states);
                end
                
                xlabel(xlab,'fontSize',12); ylabel(ylab,'fontSize',12); title(ftitle,'fontSize',13);
                
                cd('Figures')
                saveas(gca,[fname1 '.png']);
                cd('..')
            end
            
            % Make scatterplot for the correlation between percent
            % vaccination and percent change in COVID deaths
            figure
            fname2  = ['Fig-correlation-' death_types{dt} '-YEAR-' yr2 '-' strrep(ages,' ','-')]; 
            xlab   = ['Doses (' ages ', ' yr1 ') divided by population'];
            ylab   = ['% change in deaths (' yr2 ' vs ' yr1 '), ' ages];
            ftitle = ['% vaxxed vs. % change in ' death_types{dt} ' deaths: r=' num2str(out.r_per_per(i,dt),'%0.2f') ', p=' num2str(out.p_per_per(i,dt),'%0.5f')];
            
            try
                plot_fit(xvec2,yvec2,robust,'k.','r',states);
            catch
                plot_fit(xvec2,yvec2,0,'k.','r',states);
            end
            
            xlabel(xlab,'fontSize',12); ylabel(ylab,'fontSize',12); title(ftitle,'fontSize',13);
            cd('Figures')
            saveas(gca,[fname2 '.png']);
            cd('..')
        end
        
        % Estimate 2022 deaths attributed to vaccinations in 2021
        % yfit=b(1)+b(2)*dth21_log+b(3)*vax_log;
        if ~perc_reg
            % Compute VFR based on regression of raw counts
            clear numD numVax
            FF=0.1; % Sample at 10% increases

            for ii=1:length(dth1)
                Y1_log=b(1)+b(2)*log(dth1(ii))+b(3)*log(vax(ii));
                Y1=exp(Y1_log);
                Y2=Y1*exp(b(3)*log(1+FF));
                numD(ii)=Y2-Y1; numVax(ii)=FF*vax(ii);
            end
            
            rate=sum(numD)/sum(numVax);
            
            % Sample Y1_hat at median dth1 value
%             ii=find(dth1-median(dth1)==min(abs(median(dth1)-dth1)));
%             %ii=1;
%             ii=ii(1); % Just take the first sample that is near the median
%             Y1_log=b(1)+b(2)*log(dth1(ii))+b(3)*log(vax(ii));
%             Y1=exp(Y1_log);
%             Y2=Y1*exp(b(3)*log(1+FF));
%             numD=Y2-Y1; numVax=FF*vax(ii);
%             rate=numD/numVax;
            out.vfr_perc(i,dt)=rate*100;
            
            % Now compute total deaths based on estimated death rate
            out.deaths(i,dt)=rate*(sum(vax_21)+sum(vax_22));
        else
            % Compute VFR based on regression of percentages
            clear numD numVax
%             FF=0.1;
%             for ii=1:length(dth1)
%                 Y1=b(1)+b(2)*vax
%                 numD(ii)=b(2)*FF*dth1(ii);
%                 numVax(ii)=FF*vax(ii);
%             end
%             rate=sum(numD)/sum(numVax);
%             out.vfr_perc(i,dt)=rate*100;
            
%             FF=0.1; % Sample at 10% increases in percent
%             ii=1;
%             Y1=b(1)+b(2)*(vax(ii)./pop(ii));
%             Y2=b(1)+b(2)*(vax(ii)./pop(ii)+FF);
%             % convert percent change in deaths to raw deaths estimate
%             numD=(Y2-Y1)*sum(dth1); 
%             numVax=sum(vax)+FF*sum(vax);
            
            % Simple approach
            rate=b(2)*0.01*sum(dth1)/(0.01*sum(vax_21));
            out.vfr_perc(i,dt)=rate*100;
            
            % Now compute total deaths based on estimated death rate
            out.deaths(i,dt)=rate*(sum(vax_21)+sum(vax_22));
        end
        
    end
end

% Compute 95% CI 
for dt=1:length(death_types)
    out.vfr_lower_ci(1,dt)=prctile(out.vfr_perc(:,dt),2.5);
    out.vfr_upper_ci(1,dt)=prctile(out.vfr_perc(:,dt),97.5);
    out.deaths_lower_ci(1,dt)=prctile(out.deaths(:,dt),2.5);
    out.deaths_upper_ci(1,dt)=prctile(out.deaths(:,dt),97.5);
end
        

%% Here merge figures into one or more
% Load saved figures
c=hgload([fname1 '.fig']);
k=hgload([fname2 '.fig']);

% Prepare subplots
figure
h(1)=subplot(1,2,1);
h(2)=subplot(1,2,2);

% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1));
copyobj(allchild(get(k,'CurrentAxes')),h(2));

% Add legends
l(1)=legend(h(1),'LegendForFirstFigure');
l(2)=legend(h(2),'LegendForSecondFigure');

