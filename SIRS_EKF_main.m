
%Read data
Tdata = readtable('latest-ILI_incidence.csv');

%Load region data
load('regionData.mat')

%Load date file
opts = detectImportOptions('forecasting_weeks.csv');
opts.VariableTypes{1}='string';
opts.VariableTypes{3}='string';
Tdates = readtable('forecasting_weeks.csv',opts);
Tdates = Tdates(strcmp(Tdates.is_latest,'True'),:);

%Forecast date (should be Wednesday)
forecastDate = Tdates.origin_date{1};


% ============   Region-independent parameters   ============

% Rate for I -> R transition (often denoted by gamma)
pars.mu = .06;

% Rate for R -> S (loss of immunity). ILI can be caused by several
% different viruses, and even if one might become immune against some of
% them after being infected, the immunity is weaker against other viruses.
% The loss of immunity rate might not reflect on any real loss of immunity,
% but is just chosen to have a realistic time scale.
pars.phi = log(2)/60;

% Model error term to scale up the Langevin covariance
pars.CC = 4^2;

% Variance of daily change of beta
pars.Q_beta = .012^2;

pars.dnIncr = 1.05;

% ===========================================================



%Check no new countries have been introduced
CC = unique(Tdata.location);
for jc = 1:length(CC)
    if sum(strcmp(CC{jc},regionData.countryCode)) == 0
        error(['Country ' CC{jc} ' not specified'])
    end
end
cal = sort(unique(Tdata.year_week));



figure('Position',[0 0 1400 720])
Eall = zeros(2,24);
predWin = 10;
quants = [0.010 0.025 0.050 0.100 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 ... 
0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 0.975 0.990];
Tout = table;
for jc = 1:size(regionData,1)
    disp(' ')
    disp(['* * * * * * *  ' regionData.countryCode{jc} '  * * * * * * *'])
    
    %Read incidence data for the country and scale to case numbers
    Yraw = Tdata(strcmp(regionData.countryCode{jc},Tdata.location),:).value*regionData.population(jc)/1e5;
    
    
    % =============   Region-dependent parameters   =============

    %"Dark number", that is, ratio of detected and total cases. 
    pars.dn = regionData.dn(jc);
    
    % "Effective" population size (accounting for heterogeneity & observation
    % bias). Optimal for Luxembourg (with covid) was determined to be 3.2e5,
    % which is very close to half of the population. Therefore we will use
    % Nfull/2 as the population size.
    pars.N = regionData.population(jc)/2;
    
    %Coefficient for measurement error variance. Optimised for each region.
    pars.Rcoef = regionData.Rcoef(jc);
    
    
    % The dn and Rcoef parameters are set like this using data until
    % W5/2024. The coefficients 0.5 and 2 are optimised by minimising 
    % sum(Eall(1,:)).
    % pars.dn = .5*sum(Y)/length(Y)*52/pars.N; 
    % pars.Rcoef = 2*mean((Y-movmean(Y,[2 2])).^2./movmean(Y+.0001,[2 2]));
    
    % ===========================================================

    
    %Check if there is missing data and replace those by -1. Truncate the
    %data to begin from the first non-missing value
    Y = -ones(1,length(cal));
    for jw = 1:length(cal)
        ii = find(strcmp(cal{jw},Tdata(strcmp(regionData.countryCode{jc},Tdata.location),:).year_week));
        if ~isempty(ii)
            Y(jw) = Yraw(ii);
        end
    end
    Y = Y(min(find(Y>0)):length(Y));
    
    %Fill in missing data
    Y = fillData(Y);
    
    %A big outlier in the data. No effect on the projections now, but for
    %parameter tuning it should be fixed.
    if jc == 24
        Y(64) = Y(63);
    end
    
    %Edit here if you wish to go back in time and see how a projection
    %would look like
    Yfull = Y;
    Y = Y(1:end-0);  %For example Y(1:end-4) to see how the projection 4 weeks ago looked like
    
    %Run the method
    [Z,E,X,Yest,Et,dnEst] = SIRS_EKF(Y,pars);
    E(2) = sum(Et(1,:));
    E(3) = sum(Et(2,:));
    
    %Store and display results
    ii = find(Y>0);
    E0 = 0;
    for jw = ii(2:end)
        iaux = sum(ii < jw);
        E0 = E0 + abs(Y(jw)-Y(ii(iaux)))/(Y(ii(iaux))+1)^.5;
    end
    Eall(:,jc) = [E(2); E(3)];
    e10 = floor(log(min([E(1) E0]))/log(10));
    e20 = floor(log(min([E(2) E(3)]))/log(10));
    dnum(1) = E(2)/10^e20;
    dnum(2) = E(3)/10^e20;
    dnum(3) = E(1)/10^e10;
    dnum(4) = E0/10^e10;
    for jj = 1:4
        d{jj} = num2str(dnum(jj));
        if dnum(jj) < 10
            d{jj} = [' ' d{jj}];
        end
        d{jj} = [d{jj} repmat(' ',[1 10-length(d{jj})])];
    end
    
    disp('                EKF       hold')
    disp(['1-week ahead: ' d{3} d{4}])
    disp(['4-week ahead: ' d{1} d{2}])
        
    
    
    %figure;
    subplot(4,6,jc)
    hold on; grid on;
    Cin = 1e5/regionData.population(jc);
    fill([length(Yest):length(Yest)+predWin fliplr(length(Yest):length(Yest)+predWin)],Cin*[Y(end) Z(2,1:predWin) fliplr([Y(end) Z(22,1:predWin)])],[1 .87 .87],'EdgeColor','none')
    fill([length(Yest):length(Yest)+predWin fliplr(length(Yest):length(Yest)+predWin)],Cin*[Y(end) Z(7,1:predWin) fliplr([Y(end) Z(17,1:predWin)])],[1 .96 .96],'EdgeColor','none')

    plot(Cin*Y,'LineWidth',1,'Color',[0 0.4470 0.7410])
    plot(length(Y),Cin*Y(end),'ok','MarkerFaceColor','k','MarkerSize',4)
    plot(Cin*Yfull,'LineWidth',1,'Color',[0 0.4470 0.7410])
    plot(length(Yest):length(Yest)+predWin, Cin*[Y(end) Z(12,1:predWin)],'k','LineWidth',1)
    title(regionData.countryCode{jc},'FontSize',14)
    set(gca,'layer','top')
    xlim([length(Y)-10 length(Y)+predWin])
 
    
    %Create the output table for the region
    Tpred = table;
    Tpred.target = repmat({'ILI incidence'},[92 1]);
    clear('Tloc','Tor','Ttarg','Ttype')
    for jr = 1:92
        Tloc{jr,1} =  regionData.countryCode{jc};
        Tor{jr,1} = forecastDate;
        Ttarg{jr,1} = Tdates.target_end_date{Tdates.horizon==floor((jr-.001)/23)+1};
        Ttype{jr,1} = num2str(quants(1+mod(jr-1,23)));
    end
    Tpred.location = Tloc;
    Tpred.origin_date = Tor;
    Tpred.target_end_date = Ttarg;
    Tpred.output_type_id = Ttype;
    Tpred.horizon = [ones(23,1); 2*ones(23,1); 3*ones(23,1); 4*ones(23,1)];
    Tpred.output_type = repmat({'quantile'},[92 1]);
    Tpred.value = vec(Z(:,1:4)/regionData.population(jc)*1e5);
    Tpred = [Tpred; Tpred(strcmp(Tpred.output_type_id,'0.5'),:)];
    Tpred.output_type(93:96) = repmat({'median'},[4 1]);
    Tpred.output_type_id(93:96) = repmat({''},[4 1]);
    
    %Concatenate all results into one table
    Tout = [Tout; Tpred];
end
   

%sum(Eall(1,:))


% figure('Position',[300 150 650 500])
% subplot(4,1,[1 2])
% plot(X(3,:)/pars.N)
% hold on
% plot(X(1,:)/pars.N)
% grid
% legend({'Weekly new cases / N','S(t) / N'},'Location','NorthWest','FontSize',11)
% 
% subplot(4,1,3)
% plot(X(4,:).*X(1,:)/pars.N/pars.mu)
% grid
% legend({'R(t)'},'Location','NorthWest','FontSize',11)
% 
% subplot(4,1,4)
% plot(dnEst)
% grid
% legend({'Dark number'},'Location','NorthWest','FontSize',11)
% xlabel('Time (weeks)','FontSize',13)

writetable(Tout,[forecastDate '-ItaLuxColab-EpiEKF.csv'])
