function Ys=plotODEModel(options_pass,windowsize1_pass)

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

if exist('options_pass','var')==1 && isempty(options_pass)==0

    options1=options_pass;

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP, tstart1_INP, tend1_INP, printscreen1_INP]=options1();

else

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options_fit;

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
    end
end

if length(params.label)~=params.num | length(params.fixed)~=params.num | length(params.LB)~=params.num | length(params.UB)~=params.num | length(params.initial)~=params.num
    error('one or more parameter vectors do not match the number of parameters specified in <params.num>')
end


if exist('windowsize1_pass','var')==1 & isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

printscreen1=printscreen1_INP;


% <==============================================================================>
% <======================== Load epidemic data ========================================>
% <==============================================================================>

if isfile(strcat('./input/',cadfilename1,'.txt'))
    % File exists.
    data=load(strcat('./input/',cadfilename1,'.txt'));

    data=data(1:1:windowsize1,:);

else
    % File does not exist.
    data=[(0:1:windowsize1-1)' zeros(windowsize1)];
end



figure

% solve model numerically
%timevect=tstart1:1:tstart1+windowsize1-1;
timevect=data(:,1);

options=[];
IC=vars.initial;


curvess=[];
composite1=[];

SSEs=[];

for j=1:M

    param1=[];

    for i=1:params.num

        if params.fixed(i)
            param1=[param1;params.initial(i)];

        else
            param1=[param1;unifrnd(params.LB(i),params.UB(i))];
        end

    end

    if isempty(params.composite)==1
        composite1=[composite1;NaN];
    else
        composite1=[composite1;params.composite(param1')];
    end


    [~,F]=ode15s(model.fc,timevect,IC,options,param1,params.extra0);

     F=real(F);


    for i2=1:vars.num
        Ys(i2,j)={F(:,i2)};
    end


    %plot(cell2mat(Ys(1,9,:))) %M,var,time

    if vars.fit_diff==1
        fitcurve=abs([F(1,vars.fit_index);diff(F(:,vars.fit_index))]);
    else
        fitcurve=F(:,vars.fit_index);
    end

    subplot(1,2,1)
    % plot the fitting variable
    line1=plot(timevect,fitcurve,'b-')
    hold on
    subplot(1,2,2)
    line1=plot(timevect,fitcurve,'b-')
    hold on
    curvess=[curvess fitcurve];

    % compute SSE across curves
    SSEs=[SSEs;sum((fitcurve-data(:,2)).^2)];

end

[min1,index1]=min(SSEs);

subplot(1,2,1)
% plot the fitting variable
line1=plot(timevect,curvess(:,index1),'g-')
set(line1,'LineWidth',4)
hold on
subplot(1,2,2)
line1=plot(timevect,curvess(:,index1),'g-')
set(line1,'LineWidth',4)

subplot(1,2,1)
xlabel('Time')

if vars.fit_diff
    ylabel(strcat(vars.label(vars.fit_index),'''(t)',{' '}))
else
    ylabel(strcat(vars.label(vars.fit_index),'(t)',{' '}))
end

subplot(1,2,2)
xlabel('Time')
if vars.fit_diff
    ylabel(strcat(vars.label(vars.fit_index),'''(t)',{' '}))
else
    ylabel(strcat(vars.label(vars.fit_index),'(t)',{' '}))
end

%subplot(1,2,1)
%line1=plot(timevect,median(curvess,2),'k--')
%set(line1,'LineWidth',3)
%subplot(1,2,2)
%line1=plot(timevect,median(curvess,2),'k--')
%set(line1,'LineWidth',3)

% plot time series data

subplot(1,2,1)
line1=plot(data(:,1),data(:,2),'ro')
set(line1,'markersize',6,'LineWidth',2)
axis([data(1,1) data(end,1) 0 max(max(curvess))+5])

subplot(1,2,2)
line1=plot(data(:,1),data(:,2),'ro')
set(line1,'markersize',6,'LineWidth',2)
axis([data(1,1) data(end,1) 0 max(data(:,2))+5])

subplot(1,2,1)
set(gca,'FontSize', 24);
set(gcf,'color','white')
subplot(1,2,2)
set(gca,'FontSize', 24);
set(gcf,'color','white')
title('zoomed')

% plot the empirical distribution of the composite parameter

if isempty(params.composite)==0

    figure

    composite=[median(composite1) quantile(composite1,0.025) quantile(composite1,0.975)]

    cad1=strcat('\it{',params.composite_name,'}=',num2str(composite(end,1),4),' (95%CI:',num2str(composite(end,2),4),',',num2str(composite(end,3),4),')')
    hist(composite1)
    ylabel('Frequency')
    xlabel(params.composite_name)
    title(cad1)

    set(gca,'FontSize', 24);
    set(gcf,'color','white')

end


%% plot all state variables in a figure

figure

factor1=factor(vars.num);

if length(factor1)==1
    rows1=1;
    cols1=factor1;
else
    rows1=factor1(1);
    cols1=factor1(2);
end

cc1=1;
for i=1:1:vars.num

    subplot(rows1,cols1,cc1)

    plot(cell2mat(Ys(i,:,:)),'b-')
    hold on

    %for j=1:M
    %    plot(cell2mat(Ys(j,i,:)),'b-')
    %    hold on
    %end

    title(vars.label(i))
    set(gca,'FontSize', 24);
    set(gcf,'color','white')

    cc1=cc1+1;

end

for j=1:1:cols1

    subplot(rows1,cols1,rows1*cols1-cols1+j)
    xlabel('Time')
end


if 0
    % <===============================================================================================================>
    % <=========================== Generate simulated data with a given error structure ===============================>
    % <================================================================================================================>

    %dist1=0; % Normal distribution to model error structure (method1=0)
    %dist1=1; % Poisson error structure (method1=0 OR method1=1)
    %dist1=2; % Neg. binomial error structure where var = factor1*mean where
    % factor1 is empirically estimated from the time series
    % data (method1=0)
    %dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
    %dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
    %dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)

    figure

    M=1;
    dist1=1;
    factor1=1;
    d=1;

    [~,F]=ode15s(model.fc,timevect,IC,options,params.initial,params.extra0);
    if vars.fit_diff==1
        fitcurve=abs([F(1,vars.fit_index);diff(F(:,vars.fit_index))]);
    else
        fitcurve=F(:,vars.fit_index);
    end

    curves=AddErrorStructure(cumsum(fitcurve),M,dist1,factor1,d)

    plot(timevect,curves,'ro')
    hold on
    plot(timevect,fitcurve,'b-')
    curves=[(0:1:length(curves)-1)' curves];

    cad1='';
    for j=1:1:params.num
        cad1=strcat(cad1,params.label(j),'-',num2str(params.initial(j)),'-');
    end
    cad1=cell2mat(cad1);

    xlabel('Time')
    if vars.fit_diff
        ylabel(strcat('d/dt(',vars.label(vars.fit_index),'(t))'))
    else
        ylabel(strcat(vars.label(vars.fit_index),'(t)'))
    end

    set(gca,'FontSize', 24);
    set(gcf,'color','white')

    save(strcat('./input/curve-',model_INP.name,'-',cad1,'M-',num2str(M),'-dist1-',num2str(dist1),'-factor1-',num2str(factor1),'.txt'),'curves','-ascii')


end

