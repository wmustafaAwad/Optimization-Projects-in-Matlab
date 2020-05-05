

%{
clc;
clear;


verbose=1;

verbose=true;
x=[0,2];
y=[0,0,12,16];
CubicFit(x,y,~verbose)
%QuadFit(x,y,verbose)
%FibonacciOpt(region,accuracy,minimize,funct_name,verbose)
%GoldenOpt(region,accuracy,minimize,funct_name,verbose)

%}
%% Golden Section 1D-Optimization:

clc;
clear;


%Problem Parameters (Eg 6.14 , pg 356)
%f(lambda)= 200(6?4.124? +0.938?^2)(1.876? ?4.124)?1.94(3?0.97?) 
%Expected optimal step= 2.20134
syms x1 x2
f(x1,x2)=100*(x1^2-x2)^2+(1-x1)^2;
fg(x1,x2)=gradient(f); %for Cubic Only
X=[-2,-2]';
S=[0.970,0.244]';
region=[0,4]; %for Fibonacci and golden section only
normalize_s=false; %since S is already normalized from the example.

%Golden Section Parameters:
golden_params=struct();
golden_params.region=region;
golden_params.eps = 1e-2;
golden_params.minimize=1;
golden_params.verbose=false;

[X_star,f_star,num_iterations] = GoldenOpt(f,X,S,golden_params)

%% Fibonacci 1-D Optimization:
clc;
clear;


%Problem Parameters (Eg 6.14 , pg 356)
%f(lambda)= 200(6?4.124? +0.938?^2)(1.876? ?4.124)?1.94(3?0.97?) 
%Expected optimal step= 2.20134
syms x1 x2
f(x1,x2)=100*(x1^2-x2)^2+(1-x1)^2;
fg(x1,x2)=gradient(f); %for Cubic Only
X=[-2,-2]';
S=[0.970,0.244]';
region=[0,4]; %for Fibonacci and golden section only
normalize_s=false; %since S is already normalized from the example.

%Fibonacci Parameters:
fib_params=struct();
fib_params.region=region;
fib_params.eps = 1e-2;
fib_params.minimize=1;
fib_params.verbose=false;

[X_star,f_star,num_iterations] = FibonacciOpt(f,X,S,fib_params)


%% Quadratic Interpolation 
clc;
clear;
%Problem paramters: (Example from 6.10.2 Rao, example 6.9) 
%f(lambda)= lambda^2-2*lambda
%Expected optimal step= 1.00.
%Problem:
syms x1 x2
f(x1,x2)=x1-x2+2*x1^2 + 2*x1*x2 + x2^2;
X=[0,0]';
S=[-1,1]';
region=[0,4]; %For fibonacci and Golden Section Only

%QuadInterpOpt Parameters:
verbose=false;
t=0.4;
epsilon1=0.0001;
normalize_s=true;

[lambda_star,error] = SymbQuadInterpOpt(f,t,X,S,epsilon1,normalize_s,verbose)

%% Cubic Interpolation 1-D Optimization:
clc;
clear;

%Problem Parameters (Eg 6.14 , pg 356)
%f(lambda)= 200(6?4.124? +0.938?^2)(1.876? ?4.124)?1.94(3?0.97?) 
%Expected optimal step= 2.20134
syms x1 x2
f(x1,x2)=100*(x1^2-x2)^2+(1-x1)^2;
fg(x1,x2)=gradient(f); %for Cubic Only
X=[-2,-2]';
S=[0.970,0.244]';
region=[0,4]; %for Fibonacci and golden section only
normalize_s=false; %since S is already normalized from the example.


%CubicInterpOpt Parameters:
verbose=false;
t=0.25;
epsilon1=0.001;
epsilon2=0.001; %Currently not used
lmax=5;


[lambda_star,accuracyAtLambdaStar] = SymbCubicInterpOpt(f,fg,t,X, S, epsilon1, epsilon2,lmax,normalize_s,verbose)

%% Fletcher-Reeves method:
clc;
clear;

%f(x1,x2)=100*(x2-x1^2)^2+(1-x1)^2;
%NLP Problem Parametrs:



syms x1 x2 x3 x4



%Uncomment block to run
%Quadratic Example 6.9,  Rao, pg 343
f(x1,x2)=x1-x2+2*x1^2 + 2*x1*x2 + x2^2;
fg(x1,x2)=gradient(f);
X1=[0;0];
%}



%{
%Uncomment block to run
%Rosenbrock's Parabolic Valley function:
f(x1,x2)= 100*(x2-x1^2)^2 + (1-x1)^2;
X1=[-1.20000,1.00000]';
fg(x1,x2)=gradient(f);
%}  

%{
%Powell's Quartic Function

f(x1,x2,x3,x4) = (x1+ 10*x2)^2 + 5*(x3- x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4;
fg(x1,x2,x3,x4)=gradient(f);
X1=[3.0,-1.0, 0.0, 1.0]';
 %}

%{ 
%Uncomment block to run:
%Quadratic Benchmark Function (pg. 364, Rao) 
%Optimal X=[1,3]';
%Optimal F= 0;
f(x1,x2) = (x1+2*x2-7)^2 + (2*x1+x2-5)^2;
X1=[0,0]';
fg(x1,x2)=gradient(f);
%}


problem_params=struct();
problem_params.f=f;
problem_params.fg=fg;
problem_params.X1=X1;


%1D (Lambda Star) Parametrs:
lambda_star_params=struct();
lambda_star_params.verbose=false;
lambda_star_params.t=0.1;
lambda_star_params.epsilon1=0.01;
lambda_star_params.epsilon2=0.0001;
lambda_star_params.lmax=10;

%Fletcher-Reeves Parametrs:
Fletcher_Reeves_params=struct();
Fletcher_Reeves_params.eps1=1e-6;
Fletcher_Reeves_params.eps2=0.01;   
Fletcher_Reeves_params.eps3=1e-6;
Fletcher_Reeves_params.quad_interp_only=true; %set true to use quadratic interpolation only. false to use cubic primarily and quadratic when necessary.
Fletcher_Reeves_params.print_every=20; %0 to not print any function-vale. Else: Prints function value every (N) iterations.


[f_star,X_star,num_iterations] = FletcherReevesMultiOpt(problem_params,lambda_star_params,Fletcher_Reeves_params)

%% Marquardt Testing:
clc;
clear;



syms x1 x2 x3 x4

%============= Choose a Function =============%


%Example 6.9,  Rao, pg 343
f(x1,x2)=x1-x2+2*x1^2 + 2*x1*x2 + x2^2;
fg(x1,x2)=gradient(f);
H(x1,x2)=hessian(f);
%X1=[-1.2,1.0]; %1. Start with arbitrary point X1
X1=[0,0]'; %The transpose is actually important !
%}

%{
%Uncomment block to run
%Rosenbrock's Parabolic Valley function:
f(x1,x2)= 100*(x2-x1^2)^2 + (1-x1)^2;
X1=[-1.20000,1.00000]';
fg(x1,x2)=gradient(f);
H(x1,x2)=hessian(f);
%}  

%{
%Powell's Quartic Function

f(x1,x2,x3,x4) = (x1+ 10*x2)^2 + 5*(x3- x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4;
fg(x1,x2,x3,x4)=gradient(f);
H(x1,x2,x3,x4)= hessian(f);
X1=[3.0,-1.0, 0.0, 1.0]';
 %}

%== Chosen a Function ========================%

%NLP Problem Parametrs:
problem_params=struct();
problem_params.f=f;
problem_params.fg=fg;
problem_params.H=H;
problem_params.X1=X1;


%Marquardt Parametrs:
Marquardt_params=struct();
Marquardt_params.alpha=1e4;
Marquardt_params.c1=0.25;
Marquardt_params.c2=2;
Marquardt_params.eps=1e-2;



[f_star,X_star,num_iterations,num_of_inverses] = MarquardtMultiOpt(problem_params,Marquardt_params)



%% BFGS Method Testing:


clc;
clear;


syms x1 x2 x3 x4
%============= Choose a Function =============%


%Example 6.9,  Rao, pg 343
f(x1,x2)=x1-x2+2*x1^2 + 2*x1*x2 + x2^2;
fg(x1,x2)=gradient(f);
X1=[0,0]'; %The transpose is actually important !
%}


%{
%Uncomment block to run
%Rosenbrock's Parabolic Valley function:
f(x1,x2)= 100*(x2-x1^2)^2 + (1-x1)^2;
X1=[-1.20000,1.00000]';
fg(x1,x2)=gradient(f);
%}  


%{
%Powell's Quartic Function

f(x1,x2,x3,x4) = (x1+ 10*x2)^2 + 5*(x3- x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4;
fg(x1,x2,x3,x4)=gradient(f);
X1=[3.0,-1.0, 0.0, 1.0]';
%}

%== Chosen a Function ========================%


%NLP Problem Parametrs:
problem_params=struct();
problem_params.f=f;
problem_params.fg=fg;
problem_params.X1=X1;



%1D (Lambda Star) Parametrs:
lambda_star_params=struct();
lambda_star_params.verbose=false;
lambda_star_params.t=0.0001;
lambda_star_params.epsilon1=0.001;
lambda_star_params.epsilon2=0.001;
lambda_star_params.lmax=5;
lambda_star_params.normalize_s=false; 

%BFGS Parametrs:
BFGS_params=struct();
BFGS_params.B1=eye(length(X1));
BFGS_params.eps=1e-2;
BFGS_params.quad_interp_only=false;%set true to use quadratic interpolation only. false to use cubic primarily and quadratic when necessary.

[f_star,X_star,num_iterations] = BFGSMultiOpt(problem_params,lambda_star_params,BFGS_params)