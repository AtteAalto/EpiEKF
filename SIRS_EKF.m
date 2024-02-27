function [Zout,E,X,Yest,Et,dn] = SIRS_EKF(Y,pars)

% Rate for I -> R transition (often denoted by gamma)
mu = pars.mu;

% Initial rate for S -> I transition (time-dependent parameter)
beta = mu;
dn = pars.dn;

% Rate for R -> S (loss of immunity). 
phi = pars.phi;

% "Effective" population size 
N = pars.N;

% Number of initially infected
N_infected = abs(sum(Y(1:3))/2/pars.dn); 

% Variance of the size of the initially infected population
S_infected = N_infected; 

% Initial error variance of beta
S_beta = (beta/2)^2;

% Variance of daily change of beta
Q_beta = pars.Q_beta;



% State variables are: 
% X(1): S(t),  
% X(2): I(t),  
% X(3): Dummy variable keeping track of the sum of new cases over a week  
% X(4): beta(t)
Tlim = length(Y);
X = zeros(4,Tlim+1);
X(:,1) = [N-N_infected; N_infected; 0; beta];

% Initial state error covariance 
P = [S_infected -S_infected 0 0; 
    -S_infected S_infected 0 0; 
    0 0 0 0; 
    0 0 0 S_beta];

% Measurement error variance (depends on the number of infected)
R = pars.Rcoef*Y*max(1-pars.dn,.2) + pars.Rcoef*N/1e5;

% Model error term to scale up the Langevin covariance
CC = pars.CC; 

%Time series of the baseline dark number
dn = pars.dn*ones(1,Tlim);

Yest = zeros(1,Tlim); %Storage for predicted number of new cases
E = [0 0 0];
Yres = Y(1);
jweek = 0;
Et = zeros(2,Tlim);
while jweek < Tlim - .5   
    jweek = jweek + 1;

    % Number of detected cases this week depends linearly on the true number of new
    % cases this week
    C0 = [0 0 dn(jweek) 0];

    Xtemp = X(:,jweek);
    Xtemp(3) = 0;
    Phat = P;
    Phat(3,:) = 0;
    Phat(:,3) = 0;
    
    %Simulate 7 days
    for jd = 1:7
        
        %State update (Kalman filter prediction step):
        StoI = Xtemp(4)*Xtemp(2)*Xtemp(1)/N;
        ItoR = mu*Xtemp(2);
        RtoS = phi*(N-Xtemp(1)-Xtemp(2));
        
        Xtemp(1) = Xtemp(1) - StoI  + RtoS;
        Xtemp(2) = Xtemp(2) + StoI - ItoR;
        Xtemp(3) = Xtemp(3) + StoI;
        Xtemp(4) = Xtemp(4);

        %Jacobian of the dynamics function
        Jf = eye(4) + [-phi-Xtemp(4)*Xtemp(2)/N, -phi-Xtemp(4)*Xtemp(1)/N, 0, -Xtemp(2)*Xtemp(1)/N;
            Xtemp(4)*Xtemp(2)/N, -mu+Xtemp(4)*Xtemp(1)/N, 0, Xtemp(2)*Xtemp(1)/N;
            Xtemp(4)*Xtemp(2)/N, Xtemp(4)*Xtemp(1)/N, 0, Xtemp(2)*Xtemp(1)/N;
            0, 0, 0, 0];
        
        %"Process noise" covariance, assuming Langevin-type stochastics
        Q = [CC*(StoI+RtoS), -CC*StoI, -CC*StoI, 0;
            -CC*StoI, CC*(StoI + ItoR), CC*StoI, 0;
            -CC*StoI, CC*StoI, CC*StoI, 0;
            0, 0, 0, Q_beta];
        
        %Case import covariance
        Qu = (N/5000)^2*[1 -1 -1 0; -1 1 1 0; -1 1 1 0; 0 0 0 0];

        %Prediction error covariance
        Phat = Jf*Phat*Jf' + Q + Qu;
    end
    
    %Measurement covariance
    S = C0*Phat*C0'+R(jweek);
   
    %Predicted number of daily new cases
    Yest(jweek) = C0*Xtemp;
    
    %Kalman filter update step based on true and predicted number
    if Y(jweek) >= 0
        
        %Kalman gain. Make sure the gain for beta does not degenerate
        GK = Phat*C0'/S^.5;
        GK(4) = max(GK(4),pars.Q_beta^.5);
        
        %State update
        Xtemp = Xtemp + GK*(Y(jweek)-Yest(jweek))/S^.5;

        %Bounds for an acceptable solution
        Xtemp(1) = max(Xtemp(1),.3*N);
        Xtemp(2) = max(Xtemp(2),1);
        Xtemp(4) = max(Xtemp(4),.1*mu);
        
        %Store the most recent available datapoint
        Yres = Y(jweek);
        
        %Covariance update
        P = Phat - Phat*C0'*C0*Phat/S;
        
    else
        %Missing data
        P = Phat;       
    end
    X(:,jweek+1) = Xtemp;
    
    %Adapt the dark number if state estimate becomes unrealistic
    if 2*Xtemp(4)/.2 - Xtemp(1)/N/.45 > 1 
        ii = jweek-12:jweek-6;
        ii = ii(ii>0);
        dn(ii) = dn(ii).*(1+(pars.dnIncr-1)*(ii-(jweek-12))/6);

        ii = jweek-5:jweek+15;
        ii = ii(ii>0);
        ii = ii(ii < Tlim+.5);
        dn(ii) = pars.dnIncr*dn(ii);

        ii = jweek+16:jweek+26;
        ii = ii(ii < Tlim+.5);
        dn(ii) = dn(ii).*(pars.dnIncr-(pars.dnIncr-1)*(ii-(jweek+16))/10);
        
        jweek = max(jweek-12,0);
    
    %Do the 4-week ahead projection (or less) every week for parameter
    %fitting purposes
    elseif jweek < Tlim && jweek > 0
        
        Nfwd = min(Tlim-jweek,4);
        Xtemp(3) = 0;
        Ytemp = zeros(1,Nfwd);
        
        %Simulate Nfwd weeks forward
        for jj = 1:7*Nfwd
            
            StoI = Xtemp(4)*Xtemp(2)*Xtemp(1)/N;
            ItoR = mu*Xtemp(2);
            RtoS = phi*(N-Xtemp(1)-Xtemp(2));

            Xtemp(1) = Xtemp(1) - StoI  + RtoS;
            Xtemp(2) = Xtemp(2) + StoI - ItoR;
            Xtemp(3) = Xtemp(3) + StoI;
            Xtemp(4) = Xtemp(4);
            
            if mod(jj,7) == 0
                Ytemp(jj/7) = dn(jweek)*Xtemp(3);
                Xtemp(3) = 0;
            end
             
        end
        
        
        Xtemp(4) = min(Xtemp(4),3.5*mu);
        
        %Store prediction errors
        Ycomp = Y(jweek+1:jweek+Nfwd);
        
        if Et(1,jweek) == 0
            
            %Sq root
            dif1 =  (Ytemp(Ycomp>=0)).^.5-(Ycomp(Ycomp>=0)).^.5;
            dif2 = Yres^.5-(Ycomp(Ycomp>=0)).^.5;
            
            Et(1,jweek) = sum(abs(dif1));
            Et(2,jweek) = sum(abs(dif2));
            
        end
           
    end
         
end
   
E = sum(abs(Yest(Y>=0)-Y(Y>=0))./(Y(Y>=0)+1).^.5);


%% Projection for future
    
% Length (weeks) of the prediction interval
N_pred = 12;

%Quantiles to simulate
quants = [0.010 0.025 0.050 0.100 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 ... 
0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 0.975 0.990];

%Number of stochastic replicates
Nreps = 1000;

%Simulate forward in time
Z = zeros(Nreps,N_pred);
Pc = chol(P([1 2 4],[1 2 4]))';
for jscen = 1:size(Z,1) 

    rng(jscen)
    
    %Initialise estimates with the final state of the EKF
    Xtemp = X([1 2 4],end);
    Xtemp = Xtemp + Pc*randn(3,1);
    
    %Store here the predicted daily case numbers
    Ypred = zeros(1,7*N_pred);
    
    %Initialise cumulative beta-error
    for jj = 1:7*N_pred
        
        %Cumulative error in the beta-parameter accounting for fluctuations
        StoI = max(Xtemp(3)*Xtemp(2)*Xtemp(1)/N,0);
        StoI = max(StoI + pars.CC^.5*StoI^.5*randn,0);
        
        ItoR = max(mu*Xtemp(2),0);
        ItoR = max(ItoR + pars.CC^.5*ItoR^.5*randn,0);
        
        %R to S transition is not modelled in the prediction, except
        %through uncertainty
        RtoS = phi*(N-Xtemp(1)-Xtemp(2));
        RtoS = pars.CC^.5*RtoS^.5*randn;
        
        Xtemp(1) = Xtemp(1) - StoI + RtoS; 
        Xtemp(2) = Xtemp(2) + StoI - ItoR;
        Xtemp(3) = Xtemp(3) + pars.Q_beta^.5*randn;
        
        %Bounds for the solution
        Xtemp(1) = max(Xtemp(1),.3*N);
        Xtemp(2) = max(Xtemp(2),1);
        Xtemp(3) = max(Xtemp(3),.1*mu);
        Xtemp(3) = min(Xtemp(3),3.5*mu);
        
        
        Ypred(jj) = dn(end)*StoI;

    end
    
    %Store weekly data in the output matrix
    for jw = 1:N_pred
        Z(jscen,jw) = sum(Ypred((1:7) + 7*(jw-1)));
    end
    
end

%Add measurement error
Z = Z + Z.^.5.*randn(size(Z));
Z(Z<0) = 0;

%For each week, find the prediction quantiles
for jw = 1:size(Z,2)
    Z(:,jw) = sort(Z(:,jw),'ascend');
end

Zout = zeros(length(quants),size(Z,2));
for jq = 1:length(quants)
    Zout(jq,:) = sum([1; 3; 3; 1].*Z(round(Nreps*quants(jq))-1:round(Nreps*quants(jq))+2,:))/8;
end


