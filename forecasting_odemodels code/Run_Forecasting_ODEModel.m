function   [AICcs,performanceC,performanceF,forecast_model12,data1,datalatest,Ys]=Run_Forecasting_ODEModel(options_pass,tstart1_pass,tend1_pass,windowsize1_pass,forecastingperiod_pass)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% Fitting and forecasting model to epidemic data with quantified uncertainty

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

if exist('options_pass','var')==1 & isempty(options_pass)==0

    options=options_pass; %forecast horizon (number of data points ahead)

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options();

else

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options_forecast;

end

params_INP.num=length(params_INP.label); % number of model parameters

vars_INP.num=length(vars_INP.label); % number of variables comprising the ODE model


% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

DT=1;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

if method1>0
    dist1=method1;
end

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squard=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;


numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty


% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model=model_INP;
params=params_INP;
vars=vars_INP;

for j=1:params.num
    if params.initial(j)<params.LB(j) | params.initial(j)>params.UB(j)
        error('values in <params.initial> should lie within their parameter bounds defined by <params.LB> and <params.UB> ')
        return
    end
end

if length(params.label)~=params.num | length(params.fixed)~=params.num | length(params.LB)~=params.num | length(params.UB)~=params.num | length(params.initial)~=params.num
    error('one or more parameter vectors do not match the number of parameters specified in <params.num>')
end


% <==============================================================================>
% <======================== Load epidemic data ========================================>
% <==============================================================================>

data=load(strcat('./input/',cadfilename1,'.txt'));

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

if exist('forecastingperiod_pass','var')==1 & isempty(forecastingperiod_pass)==0

    forecastingperiod=forecastingperiod_pass; %forecast horizon (number of data points ahead)

else

    forecastingperiod=forecastingperiod_INP;

end

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

if exist('tstart1_pass','var')==1 & isempty(tstart1_pass)==0

    tstart1=tstart1_pass;

else
    tstart1=tstart1_INP;

end

if exist('tend1_pass','var')==1 & isempty(tend1_pass)==0

    tend1=tend1_pass;
else
    tend1=tend1_INP;

end


if exist('windowsize1_pass','var')==1 & isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

printscreen1=printscreen1_INP;

% <===========================================================================================================>
% <====== Check that the number of estimated parameters is smaller than the number of data points= ===========>
% <===========================================================================================================>

numparams=get_nparams(method1,dist1,sum(params.fixed==0),params.fixI0);

numparams
windowsize1

if numparams>=windowsize1

    error("Number of estimated parameters should be smaller than the calibration period. Consider increasing the length of the calibration period.")

end


% <==================================================================================>
% ============================ Rolling window analysis=====================================>
% <==================================================================================>

param_estims=zeros(params.num+3,3,length(tstart1:1:tend1)); % median, 95% CI: LB, UB

RMSECSS=[];
MSECSS=[];
MAECSS=[];
PICSS=[];
MISCSS=[];
RMSEFSS=[];
MSEFSS=[];
MAEFSS=[];
PIFSS=[];
MISFSS=[];

WISCSS=[];
WISFSS=[];

quantilescs=[];
quantilesfs=[];

if (tend1+windowsize1-1) > length(data(:,1))

    tend1= length(data(:,1))-windowsize1+1;
    'adjusting tend1'

end

AICcs=[];

cc1=1;
paramss=[];

composite12=[];


for i=tstart1:1:tend1  %rolling window analysis

    if printscreen1
        figure(100+i)
    end

    t_window=i:1:i+windowsize1-1;

    %calibration data
    datac=data(t_window,2);
    timevect1=DT*data(t_window,1);

    % calibration and future data
    data_all=data(i:end,2);
    timevect_all=DT*data(i:end,1);

    % first data point cannot be zero
    if datac(1)==0
        i
        warning('first data point in the time series is zero')
        
    end

    if getperformance & length(data_all)<windowsize1+forecastingperiod

        [length(data_all) windowsize1+forecastingperiod]

        error('Length of time series data is too short to evaluate the forecasting period indicated in <forecastingperiod>.')

    end


    I0=datac(1);

    data1=[data(t_window,1) datac];

    if params.fixI0
        params0=[params.initial I0 1 1];
    else
        params0=[params.initial vars.initial(vars.fit_index) 1 1];
    end

    [P_model1d,residual_model1, fitcurve_model1d, forecastcurve_model1, timevect2, initialguess,fval]=fit_model(data1,params0,numstartpoints,DT,model,params,vars,0);

    %     plot(timevect1,data1(:,2),'ko')
    %     hold on
    %     plot(timevect1, fitcurve_model1d,'r--')
    %
    %     xlabel('Time')
    %     ylabel('Cases')
    %
    %     pause


    [AICc,part1,part2,numparams]=getAICc(method1,dist1,sum(params.fixed==0),params.fixI0,fval,length(data1(:,1)))

    AICcs=[AICcs;[i AICc part1 part2 numparams]];

    if (method1==0 & dist1==1)

        factor1=1;

    elseif (method1==0 & dist1==2)  % estimate <factor1> which is the variance to mean ratio to scale the NB error structure

        binsize1=7;

        [ratios,~]=getMeanVarianceRatio(data1,binsize1,2);

        index1=find(ratios(:,1)>0);

        factor1=mean(ratios(index1,1));

        if factor1<1
            dist1=1;
            factor1=1;
        end

    elseif method1==0 & dist1==0 % normal distribution of the error structure

        var1=sum((fitcurve_model1d-datac).^2)./(length(fitcurve_model1d)-numparams); % last revised: 01 June 2022
        factor1=sqrt(var1);

    elseif method1==1
        dist1=1;
        factor1=1;

    elseif method1==3

        dist1=3;
        alpha=P_model1d(params.num+2); % VAR=mean+alpha*mean;
        factor1=alpha;

    elseif method1==4

        dist1=4;
        alpha=P_model1d(params.num+2); % VAR=mean+alpha*mean^2;
        factor1=alpha;

    elseif method1==5

        dist1=5;

        alpha=P_model1d(params.num+2); % VAR=mean+alpha*mean^d;
        d=P_model1d(params.num+3);

        factor1=alpha;

    end


    % <===========================================================================================>
    % <======= Derive parameter uncertainty of the best fitting models and save results ================================>
    % <===========================================================================================>


    fit_model1=[];
    forecast_model1=[];
    forecast_model12=[];

    Phatss_model1=zeros(M,params.num+3);

    f_model1_sims=[];

    for j=1:M


        f_model1_sim=AddErrorStructure(cumsum(fitcurve_model1d),1,dist1,factor1,d);

        f_model1_sims=[f_model1_sims f_model1_sim];


        % Fit model1 to bootstrap data
        data1=[timevect1 f_model1_sim];

        %params0=initialParams(data1(:,2),flag1);
        params0=P_model1d;

        [P_model1,residual_model1 fitcurve_model1 forecastcurve_model1 timevect2,initialguess,fval, F1,F2]=fit_model(data1,params0,2,DT,model,params,vars,forecastingperiod);

        fit_model1=[fit_model1 fitcurve_model1];

        forecast_model1=[forecast_model1 forecastcurve_model1];

        for i2=1:vars.num
            Ys(i2,j)={F2(:,i2)};
        end


        if method1==0 & dist1==0

            forecast_model12=[forecast_model12 AddErrorStructure(cumsum(forecastcurve_model1),20,dist1,factor1,0)];

        else

            forecast_model12=[forecast_model12 AddErrorStructure(cumsum(forecastcurve_model1),20,dist1,P_model1(end-1),P_model1(end))];

        end

        Phatss_model1(j,:)=P_model1;

    end %end bootstrapping loop

    data1=[timevect1 data1(:,2)];

    datalatest=data(i:end,1:2);

    % <==================================================================================================>
    % <========== Get forecast performance metrics for the model (if getperformance=1) =====================================>
    % <==================================================================================================>

    if 1

        [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,datalatest,forecast_model1,forecast_model12,forecastingperiod);

        [WISC_model1,WISFS_model1]=computeWIS(data1,datalatest,forecast_model12,forecastingperiod)

        % store metrics for calibration
        RMSECSS=[RMSECSS;[RMSECS_model1(end,end)]];
        MSECSS=[MSECSS;[MSECS_model1(end,end)]];
        MAECSS=[MAECSS;[MAECS_model1(end,end)]];
        PICSS=[PICSS;[PICS_model1(end,end)]];
        MISCSS=[MISCSS;[MISCS_model1(end,end)]];

        WISCSS=[WISCSS;[WISC_model1(end,end)]];

        % store metrics for short-term forecasts
        if forecastingperiod>0 && isempty(RMSEFS_model1)==0

            RMSEFSS=[RMSEFSS;[RMSEFS_model1(end,end)]];
            MSEFSS=[MSEFSS;[MSEFS_model1(end,end)]];
            MAEFSS=[MAEFSS;[MAEFS_model1(end,end)]];
            PIFSS=[PIFSS;[PIFS_model1(end,end)]];
            MISFSS=[MISFSS;[MISFS_model1(end,end)]];

            WISFSS=[WISFSS;[WISFS_model1(end,end)]];

        end

    end

    % <==================================================================================================>
    % <========== Compute quantiles of the calibration and forecasting periods and store ================>
    % <==================================================================================================>

    [quantilesc,quantilesf]=computeQuantiles(data1,forecast_model12,forecastingperiod);

    quantilescs=[quantilescs;quantilesc];

    quantilesfs=[quantilesfs;quantilesf];


    % <========================================================================================>
    % <================================ Parameter estimates =========================================>
    % <========================================================================================>

    % estimate median and 95% CI from distribution of parameter estimates

    for j=1:params.num

        param_estims(j,1:3,cc1) = [median(Phatss_model1(:,j)) quantile(Phatss_model1(:,j),0.025) quantile(Phatss_model1(:,j),0.975)];

    end

    param_estims(params.num+1,1:3,cc1) = [median(Phatss_model1(:,params.num+1)) quantile(Phatss_model1(:,params.num+1),0.025) quantile(Phatss_model1(:,params.num+1),0.975)];
    param_estims(params.num+2,1:3,cc1) = [median(Phatss_model1(:,params.num+2)) quantile(Phatss_model1(:,params.num+2),0.025) quantile(Phatss_model1(:,params.num+2),0.975)];
    param_estims(params.num+3,1:3,cc1) = [median(Phatss_model1(:,params.num+3)) quantile(Phatss_model1(:,params.num+3),0.025) quantile(Phatss_model1(:,params.num+3),0.975)];

    if isempty(params.composite)==0

        compositetemp=params.composite(Phatss_model1);

        composite12=[composite12;[median(compositetemp) quantile(compositetemp,0.025) quantile(compositetemp,0.975)]]

    else

        composite12=[composite12;[NaN NaN NaN]]

    end

    % Plot results

    LB1=quantile(forecast_model12',0.025)';
    LB1=(LB1>=0).*LB1;

    UB1=quantile(forecast_model12',0.975)';
    UB1=(UB1>=0).*UB1;

    median1=median(forecast_model12,2);


    % <========================================================================================>
    % <======================= Plot empirical distributions of the parameters ========================>
    % <========================================================================================>

    if printscreen1
        figure(100+i)
    end

    params1=[];
    paramslabels1=cell(1,(params.num+1)*3);

    for j=1:params.num

        if printscreen1
            subplot(2,params.num,j)
            hist(Phatss_model1(:,j))
            hold on
        end

        %line2=[param_estims(j,2,cc1)  10;param_estims(j,3,cc1)  10];
        %line1=plot(line2(:,1),line2(:,2),'r--')
        %set(line1,'LineWidth',2)

        params1=[params1 param_estims(j,1,cc1)  param_estims(j,2,cc1) param_estims(j,3,cc1) ];

        if isempty(params.label)
            xlabel(strcat('param(',num2str(j),')'))
            cad1=strcat('param(',num2str(j),')=',num2str(param_estims(j,1,cc1),2),' (95% CI:',num2str(param_estims(j,2,cc1) ,2),',',num2str(param_estims(j,3,cc1) ,2),')');

            paramslabels1(1+(j-1)*3:j*3)={strcat('param(',num2str(j),')'), strcat('param(',num2str(j),')_95%CI LB'), strcat('param(',num2str(j),')_95%CI UB')};
        else
            xlabel(params.label(j))
            cad1=strcat(cell2mat(params.label(j)),'=',num2str(param_estims(j,1,cc1),2),' (95% CI:',num2str(param_estims(j,2,cc1) ,2),',',num2str(param_estims(j,3,cc1) ,2),')');

            paramslabels1(1+(j-1)*3:j*3)={cell2mat(params.label(j)), strcat(cell2mat(params.label(j)),'_95%CI LB'), strcat(cell2mat(params.label(j)),'_95%CI UB')};
        end

        if printscreen1
            ylabel('Frequency')

            title(cad1)

            set(gca,'FontSize', 24);
            set(gcf,'color','white')
        end

    end

    params1=[params1 param_estims(j+1,1,cc1)  param_estims(j+1,2,cc1) param_estims(j+1,3,cc1) ];
    paramslabels1(1+j*3:(j+1)*3)={strcat('X0'), strcat('X0_95%CI LB'), strcat('X0_95%CI UB')};

    if method1==3 | method1==4
        params1=[params1 param_estims(j+2,1,cc1)  param_estims(j+2,2,cc1) param_estims(j+2,3,cc1) ];
        paramslabels1(1+(j+1)*3:(j+2)*3)={strcat('alpha'), strcat('alpha_95%CI LB'), strcat('alpha_95%CI UB')};

    elseif method1==5
        params1=[params1 param_estims(j+2,1,cc1)  param_estims(j+2,2,cc1) param_estims(j+2,3,cc1) ];
        paramslabels1(1+(j+1)*3:(j+2)*3)={strcat('alpha'), strcat('alpha_95%CI LB'), strcat('alpha_95%CI UB')};

        params1=[params1 param_estims(j+3,1,cc1)  param_estims(j+3,2,cc1) param_estims(j+3,3,cc1) ];
        paramslabels1(1+(j+2)*3:(j+3)*3)={strcat('d'), strcat('d_95%CI LB'), strcat('d_95%CI UB')};
    end

    % <========================================================================================>
    % <================================ Plot model fit and forecast ======================================>
    % <========================================================================================>

    if printscreen1
        subplot(2,params.num,[params.num+1:1:params.num*2])

        plot(timevect2,forecast_model12,'c')
        hold on

        % plot 95% PI

        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)

        hold on
        line1=plot(timevect2,LB1,'r--')
        set(line1,'LineWidth',2)

        line1=plot(timevect2,UB1,'r--')
        set(line1,'LineWidth',2)

        % plot model fit

        color1=gray(8);
        line1=plot(timevect1,fit_model1,'color',color1(6,:))
        set(line1,'LineWidth',1)

        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)

        % plot the data

        line1=plot(timevect_all,data_all,'bo')
        set(line1,'LineWidth',2)

        line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast_model12',0.975))*1.5];

        if forecastingperiod>0
            line1=plot(line2(:,1),line2(:,2),'k--')
            set(line1,'LineWidth',2)
        end

        axis([timevect1(1) timevect2(end) 0 max(quantile(forecast_model12',0.975))*1.5])

        xlabel('Time')

        cad2=strcat('(',caddisease,{' '},datatype,')');

        if vars.fit_diff
            ylabel(strcat(vars.label(vars.fit_index),'''(t)',{' '},cad2))
        else
            ylabel(strcat(vars.label(vars.fit_index),'(t)',{' '},cad2))
        end


        set(gca,'FontSize',24)
        set(gcf,'color','white')

        title(model.name)
    end

    % <=========================================================================================>
    % <================================ Save short-term forecast results ==================================>
    % <=========================================================================================>

    save(strcat('./output/Forecast-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'-mat')

    cc1=cc1+1;

    paramss=[paramss;params1];

end % rolling window analysis


%% plot all state variables in a figure

if vars.num>1

    figure(200)

    factor1=factor(vars.num);

    if length(factor1)==1
        rows1=1;
        cols1=factor1;
    else
        rows1=factor1(1);
        cols1=factor1(2);
    end

    cc1=1;

    for i2=1:1:vars.num

        subplot(rows1,cols1,cc1)
        %for j=1:M
        plot(quantile(cell2mat(Ys(i2,:,:))',0.5),'k-')
        hold on
        plot(quantile(cell2mat(Ys(i2,:,:))',0.025),'k--')
        plot(quantile(cell2mat(Ys(i2,:,:))',0.975),'k--')
        %end

        title(vars.label(i2))
        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        cc1=cc1+1;

    end

    for j=1:1:cols1

        subplot(rows1,cols1,rows1*cols1-cols1+j)
        xlabel('Time')
    end
end

%%


%save model parameters from tstart1 to tend1
save(strcat('./output/parameters-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'), 'param_estims','-mat')

%save calibration performance metrics from tstart1 to tend1 (AICc, MSE,
%MAE, coverage 95% PI, MIS, WIS)
save(strcat('./output/performanceCalibration-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'AICcs','RMSECSS', 'MSECSS', 'MAECSS', 'PICSS','MISCSS','WISCSS','-mat')

%save quantiles of the model fits from tstart1 to tend1
save(strcat('./output/QuantilesCalibration-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'quantilescs','-mat')

if forecastingperiod>0

    %save forecasting performance metrics from tstart1 to tend1 (AICc, MSE,
    %MAE, coverage 95% PI, MIS, WIS)
    save(strcat('./output/performanceForecasting-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'RMSEFSS', 'MSEFSS', 'MAEFSS', 'PIFSS','MISFSS','WISFSS','-mat')

    %save quantiles of the model forecasts from tstart1 to tend1
    save(strcat('./output/QuantilesForecastingPerformance-ODEModel-',cadfilename1,'-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'quantilesfs','-mat')

end

% <=============================================================================================>
% <================= Save csv file with parameters from rolling window analysis ====================================>
% <=============================================================================================>

rollparams=[(tstart1:1:tend1)' paramss];

T = array2table(rollparams);
T.Properties.VariableNames(1)={'time'};
if method1==3 | method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;
    T.Properties.VariableNames(2:(params.num+2)*3+1) = paramslabels1;
elseif method1==5   % save parameters alpha and d. VAR=mean+alpha*mean^d;
    T.Properties.VariableNames(2:(params.num+3)*3+1) = paramslabels1;
else
    T.Properties.VariableNames(2:(params.num+1)*3+1) = paramslabels1;
end

writetable(T,strcat('./output/parameters-rollingwindow-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

% <=============================================================================================>
% <================= Save csv file with composite parameter ===============================================>
% <=============================================================================================>

composite12=[(tstart1:1:tend1)' composite12];

T = array2table(composite12);
T.Properties.VariableNames(1)={'time'};
T.Properties.VariableNames(2:4) = {'composite mean','composite 95% CI LB','composite 95% CI UB'};
writetable(T,strcat('./output/parameters-composite-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

% <=====================================================================================================>
% <============================== Save file with AIC metrics ===========================================>
% <=====================================================================================================>

%[i AICc part1 part2 numparams]];

T = array2table(AICcs);
T.Properties.VariableNames(1:5) = {'time','AICc','AICc part1','AICc part2','numparams'};
writetable(T,strcat('./output/AICc-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


% <========================================================================================>
% <========================================================================================>
% <========================== Save csv file with calibration performance metrics ============================>
% <========================================================================================>
% <========================================================================================>

performanceC=[(tstart1:1:tend1)' zeros(length(MAECSS(:,1)),1)+windowsize1 MAECSS(:,1)  MSECSS(:,1) PICSS(:,1) WISCSS(:,1)];

T = array2table(performanceC);
T.Properties.VariableNames(1:6) = {'time','calibration_period','MAE','MSE','Coverage 95%PI','WIS'};
writetable(T,strcat('./output/performance-calibration-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


% <========================================================================================>
% <========================================================================================>
% <========================== Save csv file with forecasting performance metrics ============================>
% <========================================================================================>
% <========================================================================================>
 
performanceF=[];

if getperformance && forecastingperiod>0 && isempty(length(MAEFSS))==0

    performanceF=[(tstart1:1:tend1)' zeros(length(MAEFSS(:,1)),1)+forecastingperiod MAEFSS(:,1)  MSEFSS(:,1) PIFSS(:,1) WISFSS(:,1)];

    T = array2table(performanceF);
    T.Properties.VariableNames(1:6) = {'time','Horizon','MAE','MSE','Coverage 95%PI','WIS'};
    writetable(T,strcat('./output/performance-forecasting-model_name-',model.name,'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

end
