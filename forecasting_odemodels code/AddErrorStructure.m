function curves=AddErrorStructure(yi,M,dist1,factor1,d)

%yi is cumulative curve

curves=[];

for real=1:M
    
    yirData=zeros(length(yi),1);
    
    yirData(1)=yi(1);
    
    switch dist1
      
        case 0
            
            %Normal distribution
            for t=2:length(yi)
                lambda=normrnd(yi(t)-yi(t-1),factor1);
                yirData(t,1)=lambda;
            end
            
        case 1
            
            %Poisson dist
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                yirData(t,1)=poissrnd(lambda,1,1);
            end

        case 2

            % Negative binomial dist with VAR=factor*mean1;

            eps=0.001;
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                if mean1==0
                 mean1=eps;
                end

                var1=mean1*factor1;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);

            end

            % Quasi Poisson distribution

%             eps=0.001;
%             for t=2:length(yi)
%                 lambda=abs(yi(t)-yi(t-1)); % Mean of the Poisson distribution
%                 phi = factor1;    % Dispersion parameter (> 1 for overdispersion)
% 
%                 if lambda==0
%                     lambda=eps;
%                 end
% 
%                 % Step 1: Generate Gamma distributed random effects
%                 % The shape parameter of the Gamma distribution is 1/phi
%                 % The scale parameter is lambda*phi
%                 gamma_shape = 1 / phi;
%                 gamma_scale = lambda * phi;
% 
%                 gamma_samples = gamrnd(gamma_shape, gamma_scale, 1, 1);
% 
%                 % Step 2: Generate Poisson samples with the adjusted mean
%                 yirData(t,1) = poissrnd(gamma_samples);
% 
%             end


        case 3
            % Negative binomial dist with parameter VAR= MEAN + alpha*MEAN
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1+mean1*factor1;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
        case 4
            % Negative binomial dist with parameter VAR= MEAN +
            % alpha*MEAN^2
            
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1+factor1*mean1^2;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
        case 5
            % Negative binomial dist with parameter VAR= MEAN +
            % alpha*MEAN^d
            
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1+factor1*mean1^d;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
    end
    
    curves=[curves yirData];
    
end

