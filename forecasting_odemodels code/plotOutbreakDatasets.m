clear
close all

cc1=1;

for country1=[10 14 15 11 19 16 17 18]
    
    switch country1
        
     
        case 10
            
            dataname1='Zika (Antioquia)';
            
            filename1='zika-daily-onset-antioquia.txt';
            a=load(filename1);
            DT=1;
            
            tfs=20:1:60;
            
        case 11
            
            dataname1='Plague (Madagascar)';
            
            filename1='curve-plague-Madagascar-wave2.txt';
            a=load(filename1);
            DT=1;
            
            tfs=8:1:30;
            
        case 12
            
            dataname1='FMD-UK';
            
            filename1='FMD-series-UK-2001-120days.txt';
            a=load(filename1);
            DT=1;
            
            tfs=20:1:60;
            
%         case 13
%             dataname1='COVID19-Italy';
%             a=load(strcat('coronavirus-cases-country-2-05-24.txt'));
%             DT=1;
%             
%             tfs=20:1:60;
            
        case 14
            dataname1='1918 influenza (San Francisco)';
            a=load(strcat('curve-pandemic influenza-SF 1918.txt'));
            DT=1;
            tfs=20:1:42;
            
        case 15
            dataname1='2009 A/H1N1 influenza (Manitoba)';
            a=load(strcat('curve-Manitoba_wave1.txt'));
            DT=1;
            tfs=20:1:60;
      
            
        case 16
            dataname1='COVID-19 (Guangdong)';
            a=load(strcat('curve-Guangdong.txt'));
            DT=1;
            tfs=8:1:25;
            
            
        case 17
            dataname1='COVID-19 (Henan)';
            a=load(strcat('curve-Henan.txt'));
            DT=1;
            tfs=10:1:25;
            
            
%          case 18
%             dataname1='COVID19-Anhui';
%             a=load(strcat('curve-Anhui.txt'));
%             DT=1;
%             tfs=10:1:25;
            
            
       case 18
            dataname1='COVID-19 (Hunan)';
            a=load(strcat('curve-Hunan.txt'));
            DT=1;
            tfs=10:1:25;
           
        case 19
            
            dataname1='SARS (Singapore)';
            a=load('SARS-Singapore.txt');
            DT=1;
            tfs=15:1:45;

    end
    
    
    timevect1=(0:1:length(a)-1)';
    
    timevect_all=timevect1;
    
    f_ensemble_all=a(:,2);
    
    subplot(2,4,cc1)
    line1=plot(timevect1,f_ensemble_all,'b-x')
    set(line1,'LineWidth',2)
    hold on
    
    line1=[tfs(1) 0;tfs(1) max(f_ensemble_all)];
    
    line1=plot(line1(:,1),line1(:,2),'k--')
    set(line1,'LineWidth',2)

    line2=[tfs(end) 0;tfs(end) max(f_ensemble_all)];
    
    line1=plot(line2(:,1),line2(:,2),'k--')
    set(line1,'LineWidth',2)

    title(dataname1)
    axis([timevect1(1) timevect1(end) 0 max(f_ensemble_all)])
    
    if cc1>4
    xlabel('Time (days)')
    end
    
    set(gca,'FontSize',24)
    set(gcf,'color','white')
                
    cc1=cc1+1;
    
end

subplot(2,4,1)

ylabel('Cases')

    