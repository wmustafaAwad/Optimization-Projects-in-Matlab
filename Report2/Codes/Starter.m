clc;
clear;

%f(x1,x2)=100*(x2-x1^2)^2+(1-x1)^2;
%NLP Problem Parametrs:

%Example 6.9,  Rao, pg 343
syms x1 x2
f(x1,x2)=x1-x2+2*x1^2 + 2*x1*x2 + x2^2;
fg(x1,x2)=gradient(f);
X1=[0,0]'; %The transpose is actually important !

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


%%


%======================== I N P U T S =======================
%Problem Parameters:
f  = problem_params.f;                 %Multivariable Symbolic function to minimize
fg = problem_params.fg;               %Derivative of f
X1 = problem_params.X1;               %Initial point (Vector)

%Lambda Star Parametrs (Cubic Case): (Check Cubic Interpolation Function for explanation)
verbose=lambda_star_params.verbose;
t=lambda_star_params.t;
epsilon1=lambda_star_params.epsilon1;
epsilon2=lambda_star_params.epsilon2;
lmax=lambda_star_params.lmax;
normalize_s=lambda_star_params.normalize_s; 

%BFGS Parametrs:
B1=BFGS_params.B1;                 %Initial Hessian Inverse Matrix (Must be PD)
eps=BFGS_params.eps;               %To test for optimality
%====Eof: I N P U T S ========================================

%============= Prepare for the loop: ==================
Bi=B1;
Xi=X1;

exit_loop=false;
num_iterations=0;
%%%Eof Loop preparations ===============================
%%

%============================= MAIN LOOP ==================================
while(~exit_loop)
    num_iterations=num_iterations+1;
    
    %2. Calculate direction:
    grad_f_i=symToVecCalc(fg,Xi);
    Si=-Bi*grad_f_i;
    
    %3. Optimal lambda (Step size) and calculate X_(i+1)
    [lambda_star,accuracy] = SymbQuadInterpOpt(f,t,Xi,Si,epsilon1,normalize_s,verbose);
    Xi_plus1= Xi+ lambda_star.*Si;
    
    %4.Test the point Xi+1 for optimality ||GRAD(f_i+1)|| < epsilon
    grad_f_i_plus1= symToVecCalc(fg,Xi_plus1);
    
    if(sqrt(sum(grad_f_i_plus1.^2))<eps) 
        exit_major_loop=true;
        break;%Edit
    end 
    
    %5. Did not break.. So Continue updating B:
    di=lambda_star.*Si;
    gi= grad_f_i_plus1 - grad_f_i;
    diTgi= di' * gi; % di Transpose * gi (Cached since it is a common base)
    
    %Bi_new= Bi + T1 - T2 - T3  (equtaion on Rao pg 361)
    
    T1= 1+( gi' * Bi * gi   /   diTgi ) * ....
        ( di * di'   /   diTgi);
    T2= ( di * gi' * Bi )  /   diTgi ;
    T3= ( Bi * gi * di' ) /    diTgi;
    
    Bi= Bi + T1 - T2 - T3;
    
        
    
    
    
    
    %Prepare for next loop (Update) (equiv. to i=i+1):
    Xi=Xi_plus1;
    grad_f_i=grad_f_i_plus1; %Edit : Remove the top line in the loop that re-calculates this for efficiency. Becarefule
                                %not to spoil the first iteration though.
    
end

%== Eof: MAIN LOOP ========================================================


X_star =  Xi_plus1
f_star =  symToVecCalc(f,Xi_plus1)
num_iterations
