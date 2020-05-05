clc;
clear;

%Solves only problems with inequalities
%Example from  Hillier and Lieberman - Introduction to operations research pg 205

%Input problem of the form:

lin_problem=struct();


lin_problem.plot=true;


%REPLACE PROBLEM WITH YOUR OWN:
%{ 
Minimize Z = c*x,
subject to:
    Ax<=b and x >= 0

Note: b elements need not be positive .. therefore for problems with Ax>=b
    set -Ax<=-b

Note: if you want to maximize, input problem as is and set
    lin_problem.maximize=true; This just replaces c by -c.
%}

%============== REPLACE BELOW ==================================%
lin_problem.maximize=true;                  % true if you are maximizing. false if you are minimizing
lin_problem.c=[3,5];                        % 1*n vector (NOTE: note like course slides where c = n * 1)
lin_problem.A=[1,0;0,2;3,2];                % Coefficients of x in inequalities.. size= m * n
lin_problem.b=[4,12,18]';                   % Right hand side of inequalitise
%============== REPLACE ABOVE ==================================%

my_start=cputime;
[x_star, f_star, num_iterations, history] = Simplex(lin_problem);
my_end=cputime-my_start;

if(lin_problem.plot)
   scatter(0:num_iterations,history.z,'*g');
    hold on
    plot(0:num_iterations,history.z,'r');
    title(strcat('Function value vs iterations. x^*=',num2str(x_star'),', f^* = ',num2str(f_star)));
    ylabel('Function Value')
    xlabel('Number of  iterations') 
    
end


%% MATLAB Simplex:

if(lin_problem.maximize)
    c=-lin_problem.c;
else
    c=lin_problem.c;
end
A=lin_problem.A;
b=lin_problem.b;
Aeq=[];
Beq=[];
x0=0;
lb=[0,0];
ub=[inf,inf];

mat_start= cputime;
[mat_x,mat_f]=linprog(c,A,b,Aeq,Beq,lb,ub,x0);

if(lin_problem.maximize)
    mat_f=-mat_f;
end
mat_end= cputime-mat_start;


%% Print Comparison:

fprintf("\n==============================\n");
fprintf('Matlab function:\ntime= %.6f seconds\nx_star=%s\nf_star=%.6f',mat_end,mat2str(mat_x),mat_f);
fprintf("\n==============================\n");
fprintf('My function:\ntime= %.6f seconds\nx_star=%s\nf_star=%.6f\n',my_end,mat2str(x_star),f_star);

