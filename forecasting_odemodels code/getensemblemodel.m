function [curvesforecastsens,forecast1]=getensemblemodel(weights1,curvesforecastsm1,curvesforecastsm2,curvesforecastsm3)

sum(weights1)

if round(sum(weights1),3)~=1.0

    error('The sum of the weights should equal 1.0')

end

curvesforecastsens=[];

for m=1:1:length(weights1)

    switch m

        case 1

            M1=length(curvesforecastsm1(1,:));

            index1=datasample(1:M1,round(M1*weights1(m)),'Replace',false);

            if length(index1)>0
                curvesforecastsens=[curvesforecastsens curvesforecastsm1(:,index1)];
            end

        case 2

            M1=length(curvesforecastsm2(1,:));

            index1=datasample(1:M1,round(M1*weights1(m)),'Replace',false);

            if length(index1)>0
                curvesforecastsens=[curvesforecastsens curvesforecastsm2(:,index1)];
            end

        case 3

            M1=length(curvesforecastsm3(1,:));

            index1=datasample(1:M1,round(M1*weights1(m)),'Replace',false);

            if length(index1)>0
                curvesforecastsens=[curvesforecastsens curvesforecastsm3(:,index1)];
            end

    end

end

% store forecast curves
LB1=quantile(curvesforecastsens',0.025);
LB1=(LB1>=0).*LB1;

UB1=quantile(curvesforecastsens',0.975);
UB1=(UB1>=0).*UB1;

forecast1=[median(curvesforecastsens,2) LB1' UB1'];

