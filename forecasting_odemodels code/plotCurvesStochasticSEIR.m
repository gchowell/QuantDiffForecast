clear
close all

%curves=load('curves-SEIR-2.6-1.3-tc-30-N-10000.txt');

curves=load('curves-SEIR-2-1-0.5-tc-30-tc-40-N-100000.txt');


cc1=1;

%for country1=[10:12 14:19]
for i=[1 2 4 5 ]
    
    dataname1=strcat('SEIR-2-1-0.5-tc-30-tc-40-N-100000-curve-',num2str(i));

    
    a=[(0:1:length(curves)-1)' curves(:,i)];
    DT=1;
    
    %tfs=15:1:35;
         
    tfs=20:1:40;
    
    
    
    timevect1=(0:1:length(a)-1)';
    
    timevect_all=timevect1;
    
    f_ensemble_all=a(:,2);
    
    subplot(2,2,cc1)
    
    plot(timevect1,f_ensemble_all,'bo')
    hold on
    
    
    
    line1=[15 0;15 120];
    
    line2=plot(line1(:,1),line1(:,2),'k--')
    set(line2,'LineWidth',2)
    
    line1=[40 0;40 120];
    
    line2=plot(line1(:,1),line1(:,2),'k--')
    set(line2,'LineWidth',2)
    
    axis([0 60 0 120])
    
    xlabel('Time (days)')
    ylabel('Cases')
    
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    
    cc1=cc1+1;
    
end

