
%Code To Compare the Multi-Variable Optimization Techniques. 

%All the techniques share a similar interface, so you can adjust the problem as you may and watch the results without editing 
    %any method-related paramters.

clc;
clear;
%% Problem Definition:

syms x1 x2 x3 x4

%====================================================================================%
%$$$$$$$$$$$          EDIT PART BELOW TO TEST        $$$$$$$$$$$$%
%Example 6.9,  Rao, pg 343
%f(x1,x2)=x1-x2+2*x1^2 + 2*x1*x2 + x2^2;     %Function to optimize
%fg(x1,x2)=gradient(f);                      %Gradient of the function
%H(x1,x2)=hessian(f);                        %Hessian of the Function. Useful for Marquardt Only.
%X=[0,0]';                                   %Initial X. The transpose is actually important !


%Rosenbrock's Parabolic Valley function:
f(x1,x2)= 100*(x2-x1^2)^2 + (1-x1)^2;       
X=[-1.20000,1.00000]';
fg(x1,x2)=gradient(f);
H(x1,x2)=hessian(f);
%}


%{
%Powell's Quartic Function

f(x1,x2,x3,x4) = (x1+ 10*x2)^2 + 5*(x3- x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4;
fg(x1,x2,x3,x4)=gradient(f);
H(x1,x2,x3,x4)= hessian(f);
X=[3.0,-1.0, 0.0, 1.0]';
 %}
%$$$$$$$$$$$         EDIT PART ABOVE TO TEST        $$$$$$$$$$$$$%
%====================================================================================%

%Assign NLP Problem Parametrs:
problem_params=struct();
problem_params.f=f;
problem_params.fg=fg;
problem_params.H=H;
problem_params.X1=X;


%% Problem-irrelevant parameters for methods:
verbose=false;

%1D (Lambda Star) Parametrs:
lambda_star_params=struct();
lambda_star_params.verbose=false;
lambda_star_params.t=0.1;
lambda_star_params.epsilon1=0.001;
lambda_star_params.epsilon2=0.001;
lambda_star_params.lmax=5;
lambda_star_params.normalize_s=false; 

%Fletcher-Reeves Parametrs:
Fletcher_Reeves_params=struct();
Fletcher_Reeves_params.eps1=1e-6;
Fletcher_Reeves_params.eps2=0.01;   
Fletcher_Reeves_params.eps3=1e-6;
Fletcher_Reeves_params.quad_interp_only=true; %set true to use quadratic interpolation only. false to use cubic primarily and quadratic when necessary.
Fletcher_Reeves_params.print_every=0; %0 to not print any function-vale. Else: Prints function value every (N) iterations.


%Marquardt Parametrs:
Marquardt_params=struct();
Marquardt_params.alpha=1e4;
Marquardt_params.c1=0.25;
Marquardt_params.c2=2;
Marquardt_params.eps=1e-2;

%BFGS Parametrs:
BFGS_params=struct();
BFGS_params.B1=eye(length(X));
BFGS_params.eps=1e-2;
BFGS_params.quad_interp_only=true;%set true to use quadratic interpolation only. false to use cubic primarily and quadratic when necessary.



%% Prepare:

%All functions output: [f_star,X_star,num_iterations]
display_names=["Fletcher-Reeves method",...
              "Marquardt Method", ...
              "BFGS (Quasi-Newton) Method",...
              ];
          
method_functions=["FletcherReevesMultiOpt", ... 
                   "MarquardtMultiOpt", ... 
                   "BFGSMultiOpt",...
                   ];
               
inputs=[        "(problem_params,lambda_star_params,Fletcher_Reeves_params)",...
                "(problem_params,Marquardt_params)" ,...
                "(problem_params,lambda_star_params,BFGS_params)",...
                 ];
             
X_star=0;
f_star=0;
num_iterations=0;



%% Compare:

for i=1:length(display_names)
    time=0;
    start_t=cputime; 
    [f_star,X_star,num_iterations]=eval(strcat(method_functions(i),inputs(i)));
    time=cputime-start_t;
    
    fprintf("=============\n%s:\n-Optimal X:\t\t\t\t\t%s\n-num_iterations:\t\t\t%d\n-Optimal Function_Value:\t%.10f\n-cputime:\t\t\t\t\t%.10f seconds\n=============\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",display_names(i),mat2str(X_star),num_iterations,f_star,time);
    
    X_star=0;
    f_star=0;
    num_iterations=0;
    
end
