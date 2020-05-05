
%Code To Compare the Single Variables Optimization Techniques. If a function produces an error, it is skipped. You can see that from the number of iterations, which will be zero.

%All the techniques share a similar interface, so you can adjust the problem as you may and watch the results. 

%Please note that the cubic-interpolation method fails (produces an error) in case of quadratic function of lambda.
clc;
clear;
%% Problem Definition:
%Problem Parameters (Eg 6.14 , pg 356)
%f(lambda)= 200(6?4.124? +0.938?^2)(1.876? ?4.124)?1.94(3?0.97?)  %?=lambda
%Expected optimal step= 2.20134

%====================================================================================%
%$$$$$$$$$$$          EDIT PART BELOW TO TEST        $$$$$$$$$$$$%

syms x1 x2                                          %Make sure to define all used symbolic variables here

f(x1,x2)=100*(x1^2-x2)^2+(1-x1)^2;                  %Function to optimize
fg(x1,x2)=gradient(f);                              %Useful for Cubic-Interpolation Only

X=[-2,-2]';                                         %Starting point (lambda=0 point)
S=[0.970,0.244]';                                   %Descent Direction for unimodal function
%$$$$$$$$$$$         EDIT PART ABOVE TO TEST        $$$$$$$$$$$$$%
%====================================================================================%

region=[0,4];                                       %Search region. Useful for Fibonacci and golden section only.
normalize_s=false;                                  %since S is already normalized from the example.
t=0.4;                                              %Step size in point-search. Useful in Quad and Cubic only.



%% Problem-irrelevant parameters for methods:
verbose=false;

%Golden Section Parameters:
golden_params=struct();
golden_params.region=region;
golden_params.eps = 1e-2;
golden_params.minimize=1;
golden_params.verbose=false;

%Fibonacci Parameters:
fib_params=struct();
fib_params.region=region;
fib_params.eps = 1e-2;
fib_params.minimize=1;
fib_params.verbose=false;

%Cubic and Quadratic Parameters:
epsilon1=0.001;
epsilon2=0.001; %Currently not used
lmax=5;



%% Prepare:

%All functions output: [lambda_star, __ ,  num_iterations/num_refits]
display_names=["Fibonacci method",...
              "Golden Section Method", ...
              "Quadratic-Interpolation Method",...
              "Cubic-Interpolation Method"];
          
method_functions=["FibonacciOpt", ... 
                   "GoldenOpt", ... 
                   "SymbQuadInterpOpt",...
                   "SymbCubicInterpOpt"...
                   ];
               
inputs=[        "(f,X,S,fib_params)",...
                "(f,X,S,golden_params)" ,...
                "(f,t,X,S,epsilon1,normalize_s,verbose)",...
                "(f,fg,t,X, S, epsilon1, epsilon2,lmax,normalize_s,verbose)"....
                 ];
             
X_stars=zeros(length(display_names));
f_stars=zeros(length(display_names));
nums_iterations=zeros(length(display_names));


%% Compare:

for i=1:length(display_names)
    time=0;
    start_t=cputime; 
    f_star=0;
    [X_stars(i),discarded_var,nums_iterations(i)]=eval(strcat(method_functions(i),inputs(i)));
    f_stars(i)=symToVecCalc(f,X+X_stars(i).*S);
    time=cputime-start_t;
    
    fprintf("=============\n%s:\n-Lambda*:\t\t\t%.10f\n-num_iterations:\t%d\n-Function_Value:\t%.10f\n-cputime:\t\t\t%.10f seconds\n=============\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",display_names(i),X_stars(i),nums_iterations(i),f_stars(i),time);
    
    
end
