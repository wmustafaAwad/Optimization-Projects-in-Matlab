
clc;
clear;

%Assume all are equalities
%========================== Problem Paramteres ========================


%{
%Uncomment block to run
%Problem 1:
problem=struct();
problem.maximize=true;
problem.c=[1.1;1;0];
problem.A=[1,1,1];
problem.b=[6];

%}


%Problem 2:
problem=struct();
problem.maximize=true;          % true if maximization, false if minimization
problem.c=[30;20;0;0];          % Weights of x's in objective function
problem.A=[2,1,1,0;1,3,0,1];    % m * n matrix. Weights of x's in m constraint equality equations
problem.b=[8;8];                % m*1 vector (column vector)

%}


problem.gamma= 1e-2;                % Used in long-step method (0,1)
problem.sigma=0.5;                  % [0:1)(Centering parameter)
problem.StopMu=1e-3;                % Stop when mu is <= StopMu
problem.plot=true;                 % If true, history is plotted
problem.plot_stop_for=5;            % Stop plot for 5 seconds
problem.auto_init=true;             % Heuristic initialization? (otherwise: x=ones, s= ones, lambda=zeros)
problem.auto_init_diversify=false;  % Diversify during automatic initialization ? (Not done in slides example)
problem.save_figures=true;          % True to save resulting figures

if(problem.auto_init)
    [problem.init_x,problem.init_s ,problem.init_lambda] = IPinitParams(problem);
end
%== Eof: Problem Paramteres ============================================



model_names=["Primal-Dual Path-Following","Long-step Path-Following","Predictor-Corrector (Mehrotra)"];

model_functs=["IpPrimDualPathFollowing","IpLongStepPathFollowing","IpPredictorCorrector"];




for i=1:length(model_names)
    close all;
    [x_star,f_star,history]=eval(strcat(model_functs(i),'(problem);'));
    iters=history.iters-1; %One step is zero step
    
    fprintf("\n================\nFor %s:\n.:x_star=%s\n.:f_star=%.6f\n.:n_iterations=%d",model_names(i),mat2str(x_star'),f_star,iters);
    
    if(problem.plot)
        DrawIPHistory(history,problem,model_names(i));
    end

    fprintf("\nPress Enter to Continue .. %d methods left to try\n",length(model_names)-i);
    pause;
    %Delete last line:
    fprintf(repmat('\b', 1, length('\Press Enter to Continue .. 9 methods left to try\')));
end
fprintf('\n');



