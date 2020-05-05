function [x_init,s_init,lambda_init] = IPinitParams(problem)


    AATinv=inv((problem.A*problem.A'));
    problem.init_x=     problem.A'*AATinv*problem.b;            %Initial x ..
    problem.init_lambda= AATinv*problem.A*problem.c;            %Initial lambda ...
    problem.init_s= problem.c-problem.A'*problem.init_lambda;   %Initial s ...=

    e=ones(size(problem.A,2),1);
    problem.init_x = problem.init_x + max(-1.5*min(problem.init_x),0)*e;
    problem.init_s = problem.init_s + max(-1.5*min(problem.init_s),0)*e;

    sigma_x=0.5 * (problem.init_x' * problem.init_s) / (e' * problem.init_s);
    sigma_s=0.5 * (problem.init_x' * problem.init_s) / (e' * problem.init_x);

    try
        if(problem.auto_init_diversify)
            problem.init_x = problem.init_x + sigma_x*e;
            problem.init_s = problem.init_s + sigma_s*e;
        end
    catch
        nothing=0;
    end
    
    
    x_init      = problem.init_x;
    s_init      = problem.init_s;
    lambda_init = problem.init_lambda;









end