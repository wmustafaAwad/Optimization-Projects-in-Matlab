function [f_star,X_star,num_iterations,num_of_inverses] = MarquardtMultiOpt(problem_params,Marquardt_params)
    %Mutivariable optimization via the Marquardt Newton Method

    %======================== I N P U T S =======================
    %Problem Parameters:
    f  = problem_params.f;                %Multivariable Symbolic function to minimize
    fg = problem_params.fg;               %Derivative of f
    H  = problem_params.H;                %Hessian Matrix of f
    X1 = problem_params.X1;               %Initial point (Vector)


    %Marquardt Parametrs:
    alpha=Marquardt_params.alpha;           %To maintain Positive-Definiiteness
    c1=Marquardt_params.c1;                 %To decerease alpha
    c2=Marquardt_params.c2;                 %To increase alpha
    eps=Marquardt_params.eps;               %To test for optimality
    %====Eof: I N P U T S ========================================


    %============= Prepare for the loop: ==================
    exit_major_loop=false;
    exit_minor_loop=false;
    Xi=X1;
    fi=symToVecCalc(f,Xi);
    num_iterations=0;
    num_of_inverses=0;
    %%%Eof Loop preparations ===============================
    %%
    %2.  gradient f (xi):
    grad_f_xi=symToVecCalc(fg,Xi);

    %============================= MAIN LOOP ==================================
    while(~exit_major_loop)
        num_iterations=num_iterations+1;
        %3. Optimal?
        if(sqrt(sum(grad_f_xi.^2))<eps) %If so, it's optimal 
            exit_major_loop=true;
            break;%Edit
        else %4. Not Optimal.. go on
            Si= -1 * double(symToVecCalc(H,Xi)+alpha.*eye(length(Xi)))\symToVecCalc(fg,Xi); num_of_inverses=num_of_inverses+1;

            Xi_test=Xi+Si;
            fi_plus1=symToVecCalc(f,Xi);

            while(~exit_minor_loop)
                if(fi_plus1<fi) %5. Compare and Move to next Step !
                    alpha=c1*alpha; %6. decrease alpha
                    exit_minor_loop=true; 
                else %Hmmm.. need positive definit step.. let's make alpha larger
                    alpha=c2*alpha; %7. Increase alpha (Bring on the positive-definitiess)
                    Si= -1 * double(symToVecCalc(H,Xi)+alpha.*eye(length(Xi)))\symToVecCalc(fg,Xi); num_of_inverses=num_of_inverses+1;
                    Xi_test=Xi+Si;
                    fi_plus1=symToVecCalc(f,Xi_test);
                end
            end

            Xi=Xi_test;
            fi=fi_plus1;
            grad_f_xi=symToVecCalc(fg,Xi);
            exit_minor_loop=false;

        end
    end

    %== Eof: MAIN LOOP ========================================================

    X_star=Xi;
    f_star=fi;
    %num_iterations
    %num_of_inverses




end
