function numparams=get_nparams(method1,dist1,params_num,fixI0)

numparams=params_num;

numparams=numparams+fixI0;


if method1==0 & dist1==0 % Normal distribution -- one parameter for variance

    numparams=numparams+1;

elseif method1==0 & dist1==2 % Quasi Poisson distribution -- one parameter for the dispersion parameter

    numparams=numparams+1;

elseif method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)

    numparams=numparams+1;

elseif method1==5

    numparams=numparams+2;  %Neg. Binomial requires 2 more parameters (alpha,d)

end


