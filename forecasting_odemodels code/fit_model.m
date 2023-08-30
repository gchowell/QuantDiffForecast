% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [P residual fitcurve forecastcurve timevect2,initialguess,fval,F1,F2]=fit_model(data1,params0,numstartpoints,DT,modelX,paramsX,varsX,forecastingperiod)

global model params vars method1 timevect ydata

model=modelX;
params=paramsX;
vars=varsX;

timevect=data1(:,1)*DT;

%timevect=(data1(:,1))*DT;

I0=data1(1,2); % initial condition

z=params0;

for i=1:params.num

    if params.fixed(i)

        params.LB(i)=params0(i);

        params.UB(i)=params0(i);

    end

end

switch method1

    case 0
        LBe=[0 0];
        UBe=[0 0];
    case 1
        LBe=[0 0];
        UBe=[0 0];
    case 3
        LBe=[10^-8 1];
        UBe=[10^5 1];
    case 4
        LBe=[10^-8 1];
        UBe=[10^5 1];
    case 5
        LBe=[10^-8 0.6]; %d>=1
        UBe=[10^5 10^3];

        %LBe=[10^-8 0.5];
        %UBe=[10^5 0.5];

end

if params.fixI0==1

    LB=[params.LB I0 LBe];
    UB=[params.UB I0 UBe];

else

    LB=[params.LB 0 LBe];
    UB=[params.UB sum(abs(data1(:,2))) UBe];

end

% 
% if 0 % USE LSQCURVEFIT (Non-linear least squares)
% 
%     options=optimset('tolfun',10^-5,'TolX',10^-5,'MaxFunEvals',3200,'MaxIter',3200, 'algorithm','trust-region-reflective');
% 
%     [P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowth1,z,timevect,data1(:,2),LB,UB,options);
% 
%     f=@plotModifiedLogisticGrowth1;
% 
%     problem = createOptimProblem('lsqcurvefit','x0',z,'objective',f,'lb',LB,'ub',UB,'xdata',timevect,'ydata',data1(:,2),'options',options);
% 
%     %ms = MultiStart('PlotFcns',@gsplotbestf,'Display','final');
% 
%     ms = MultiStart('Display','final');
% 
%     ms = MultiStart(ms,'StartPointsToRun','bounds')
% 
%     [P,errormulti] = run(ms,problem,20)
% 
%     z=P;
% 
%     [P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowth1,z,timevect,data1(:,2),LB,UB,options);
% 
% 
% end

%A=[];      % We are using fmincon, but using none of the constraint options
%b=[];
%Aeq=[];
%beq=[];
%nonlcon=[];

%options=optimset('tolfun',10^-5,'TolX',10^-5,'MaxFunEvals',3200,'MaxIter',3200, 'algorithm','interior-point');

%[P, fval, exitflag]=fmincon(@plotModifiedLogisticGrowthMethods1,z,A,b,Aeq,beq,LB,UB,nonlcon,options);

%method1=1; %LSQ=0, MLE (Poisson)=1, Pearson chi-squared=2. MLE(neg binomial)=3

ydata=data1(:,2);

options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);

%options=optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',10000,'MaxIter',10000);

f=@parameterSearchODE;

problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);

%ms = MultiStart('PlotFcns',@gsplotbestf);
%ms = MultiStart('Display','final');
ms = MultiStart('Display','off');

%pts = z;
tpoints = CustomStartPointSet(z);

flagg=-1;

while flagg<0

    initialguess=[];

    if numstartpoints>0
        rpoints = RandomStartPointSet('NumStartPoints',numstartpoints); % start with a few random starting sets in addition to the guess supplied by the user (z)

        allpts = {rpoints,tpoints};
        initialguess=list(rpoints,problem);

    else
        allpts = {tpoints};

    end

    initialguess=[initialguess;z];

    %z
    %list(tpoints)

    %ms = MultiStart(ms,'StartPointsToRun','bounds')
    %[xmin,fmin,flag,outpt,allmins] = run(ms,problem,allpts);

    [P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);

end

% ydata
% initialguess
% P
% pause

% P is the vector with the estimated parameters

options = [];

IC=vars.initial;

if params.fixI0==1
    IC(vars.fit_index)=I0;
else
    IC(vars.fit_index)=P(params.num+1);
end

[~,F]=ode15s(model.fc,timevect,IC,options,P,params.extra0);
F1=F;

if vars.fit_diff==1
    fitcurve=abs([F(1,vars.fit_index);diff(F(:,vars.fit_index))]);
else
    fitcurve=F(:,vars.fit_index);
end

residual=fitcurve-ydata;

%fitcurve=residual+data1(:,2);

if forecastingperiod<1

    forecastcurve=residual+data1(:,2);
    timevect2=timevect;
    F2=F1;

else
    timevect2=(data1(1,1):data1(end,1)+forecastingperiod)*DT;

    [~,F]=ode15s(model.fc,timevect2,IC,options,P,params.extra0);
    F2=F;

    if vars.fit_diff==1
        forecastcurve=abs([F(1,vars.fit_index);diff(F(:,vars.fit_index))]);
    else
        forecastcurve=F(:,vars.fit_index);
    end

end
