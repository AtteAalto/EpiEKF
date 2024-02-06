function [Z,E,X,Yest,Et] = SIRS_EKF(Y,pars)

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
R = pars.Rcoef*Y*(1-pars.dn) + pars.Rcoef*N/1e5;

% Model error term to scale up the Langevin covariance
CC = pars.CC; 

% Number of detected cases this week depends linearly on the true number of new
% cases this week
C0 = [0 0 pars.dn 0];

Yest = zeros(1,Tlim); %Storage for predicted number of new cases
E = [0 0 0];
Yres = Y(1);
for jweek = 1:Tlim
    
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
        Xtemp(4) = min(Xtemp(4),3.5*mu);
        
        %Store the most recent available datapoint
        Yres = Y(jweek);
        
        %Covariance update
        P = Phat - Phat*C0'*C0*Phat/S;
        
    else
        %Missing data
        P = Phat;       
    end

    X(:,jweek+1) = Xtemp;
    
    %Do the 4-week ahead projection (or less) every week for parameter
    %fitting purposes
    if jweek < Tlim
        
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
                Ytemp(jj/7) = dn*Xtemp(3);
                Xtemp(3) = 0;
            end
             
        end
        
        %Store prediction errors
        Ycomp = Y(jweek+1:jweek+Nfwd);
        Et(1,jweek) = sum((Ytemp(Ycomp>0)-Ycomp(Ycomp>0)).^2);
        Et(2,jweek) = sum((Yres-Ycomp(Ycomp>0)).^2);
        E(2) = E(2) + sum((Ytemp(Ycomp>0)-Ycomp(Ycomp>0)).^2);
        E(3) = E(3) + sum((Yres-Ycomp(Ycomp>0)).^2);
             
    end
        
end
   
E(1) = sum((Yest(Y>0)-Y(Y>0)).^2);


%% Projection for future
    
% Length (weeks) of the prediction interval
N_pred = 10;

%Quantiles to simulate
quants = [0.010 0.025 0.050 0.100 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 ... 
0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 0.975 0.990];
quants = [quants, (quants(1:end-1)+quants(2:end))/2, .005 .995];
quants = sort(quants,'ascend');

%From quantiles to beta-coefficients
coefs = 2^.5*erfinv(2*quants-1);

%Simulate forward in time
Z = zeros(length(quants),N_pred);
for jscen = 1:size(Z,1) 

    %Initialise estimates with the final state of the EKF
    Xtemp = X([1 2 4],end);
    
    %Store here the predicted daily case numbers
    Ypred = zeros(1,7*N_pred);
    
    %Initialise cumulative beta-error
    cumerr = P(4,4);
    for jj = 1:7*N_pred
        
        %Cumulative error in the beta-parameter accounting for fluctuations
        cumerr = cumerr + pars.Q_beta;
        
        StoI = max(Xtemp(3)+coefs(jscen)*cumerr^.5,0)*Xtemp(2)*Xtemp(1)/N;
        Xtemp(1) = Xtemp(1) - StoI; 
        Xtemp(2) = Xtemp(2) + StoI - mu*Xtemp(2);
        Ypred(jj) = dn*StoI;

    end
    
    %Store weekly data in the output matrix
    for jw = 1:N_pred
        Z(jscen,jw) = sum(Ypred((1:7) + 7*(jw-1)));
    end
    
end

Ypred = Z((size(Z,1)+1)/2,:);
for jw = 1:size(Z,2)
    Z(:,jw) = sort(Z(:,jw),'ascend');
end
Z = Z(2:2:end,:);

%Add measurement error
for jscen = 1:size(Z,1)
    Z(jscen,:) = Z(jscen,:) + coefs(jscen)*Ypred.^.5;
end

%No negative projections
Z(Z<0) = 0;




  