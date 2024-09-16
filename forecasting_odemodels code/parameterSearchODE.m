function objfunction=parameterSearchODE(z)

global model params vars method1 timevect ydata

numFitIndices = length(vars.fit_index);  % Cache the number of fitting indices for clarity and efficiency

I0=z(params.num+1:params.num+numFitIndices);
alpha=z(params.num+numFitIndices+1);
d=z(params.num+numFitIndices+2);


% Initialize initial conditions from predefined variables
IC = vars.initial;

% Set specific initial conditions for indices specified in vars.fit_index
IC(vars.fit_index) = I0;

% Solve the ODE using ode15s; here 'model.fc' is the model function, 'timevect' is the time vector
[~, x] = ode15s(model.fc, timevect, IC, [], z, params.extra0);

% Initialize yfit array to store fit results
yfit = zeros(length(ydata), 1);

% Initialize an index to keep track of the end position in yfit
currentEnd = 0;

% Loop through each fit index to process the ODE solution
for j = 1:length(vars.fit_index)
    % Check if differentiation is needed based on vars.fit_diff
    if vars.fit_diff(j) == 1
        % Calculate the absolute derivative of the ODE solution for this variable
        fitcurve = abs([x(1, vars.fit_index(j)); diff(x(:, vars.fit_index(j)))]);
    else
        % Use the ODE solution directly
        fitcurve = x(:, vars.fit_index(j));
    end

    % Append the fitcurve to the yfit array, updating currentEnd to the new position
    yfit(currentEnd + 1 : currentEnd + length(fitcurve)) = fitcurve;
    currentEnd = currentEnd + length(fitcurve);
end


eps=0.001;

%%MLE expression
%This is the negative log likelihood, name is legacy from least squares code
%Note that a term that is not a function of the params has been excluded so to get the actual
%negative log-likliehood value you would add: sum(log(factorial(sum(casedata,2))))

if sum(yfit)==0
    objfunction=10^10;%inf;
else
    %    z
    yfit(yfit==0)=eps; %set zeros to eps to allow calculation below.  Shouldn't affect solution, just keep algorithm going.


    switch method1

        case 0  %Least squares


            %SSE=sum((ydata-yfit).^2);

            %objfunction= -1*((-length(ydata)/2)*log(2*pi)-(length(ydata)/2)*log(SSE/length(ydata))-length(ydata)/2);


            objfunction=sum((ydata-yfit).^2);



        case 1 % MLE for Poisson distribution (negative log-likelihood)


            %             sum1=0;
            %             for i=1:length(ydata)
            %
            %                 sum1=sum1+ydata(i)*log(yfit(i))-sum(log(2:1:ydata(i)))-yfit(i);
            %
            %             end
            %
            %             objfunction=-sum1;

            objfunction=-sum(ydata.*log(yfit)-yfit);


        case 3  % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean;


            sum1=0;

            for i=1:length(ydata)
                for j=0:(ydata(i)-1)

                    sum1=sum1+log(j+(1/alpha)*yfit(i));

                end

                %sum1=sum1+ydata(i)*log(alpha)-(ydata(i)+(1/alpha)*yfit(i))*log(1+alpha)-sum(log(2:1:ydata(i)));
                sum1=sum1+ydata(i)*log(alpha)-(ydata(i)+(1/alpha)*yfit(i))*log(1+alpha);

            end

            objfunction=-sum1;

        case 4
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^2;

            sum1=0;

            for i=1:length(ydata)
                for j=0:(ydata(i)-1)

                    sum1=sum1+log(j+(1/alpha));

                end

                %sum1=sum1+ydata(i)*log(alpha*yfit(i))-(ydata(i)+(1/alpha))*log(1+alpha*yfit(i))-sum(log(2:1:ydata(i)));
                sum1=sum1+ydata(i)*log(alpha*yfit(i))-(ydata(i)+(1/alpha))*log(1+alpha*yfit(i));

            end

            objfunction=-sum1;

        case 5
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^d;

            sum1=0;

            for i=1:length(ydata)
                for j=0:(ydata(i)-1)

                    sum1=sum1+log(j+(1/alpha)*yfit(i).^(2-d));

                end

                %sum1=sum1+ydata(i)*log(alpha*(yfit(i).^(d-1)))-(ydata(i)+(1/alpha)*yfit(i).^(2-d))*log(1+alpha*(yfit(i).^(d-1)))-sum(log(2:1:ydata(i)));
                sum1=sum1+ydata(i)*log(alpha*(yfit(i).^(d-2)).*yfit(i))-(ydata(i)+(1/alpha)*yfit(i).^(2-d))*log(1+alpha*(yfit(i).^(d-2)).*yfit(i));

            end

            objfunction=-sum1;

        case 6

           objfunction=sum(abs(ydata-yfit));


    end

end
