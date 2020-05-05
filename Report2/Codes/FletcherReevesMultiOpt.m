function [f_star,X_star,num_iterations] = FletcherReevesMultiOpt(problem_params,lambda_star_params,Fletcher_Reeves_params)
    %Conjugate Gradient Multivariable Optimization via FLetcher-Reeves
    %method (Rao, Chapter 6)
    %======================== I N P U T S =======================
    %Problem Parameters:
    f=problem_params.f;                 %Multivariable Symbolic function to minimize
    fg=problem_params.fg;               %Derivative of f
    X1=problem_params.X1;               %Initial point (Vector)

    %Lambda Star Parametrs (Cubic Case): (Check Cubic Interpolation Function for explanation)
    verbose=lambda_star_params.verbose;
    t=lambda_star_params.t;
    epsilon1=lambda_star_params.epsilon1;
    epsilon2=lambda_star_params.epsilon2;
    lmax=lambda_star_params.lmax;

    %Fletcher-Reeves Parametrs:
    eps1=Fletcher_Reeves_params.eps1;   %For function  convergence condition
    eps2=Fletcher_Reeves_params.eps2;   %For derivative  convergence condition
    eps3=Fletcher_Reeves_params.eps3;   %For X convergence condition
    %====Eof: I N P U T S ========================================


    %First Step:
    % 2. Set the first direction S1 to -gradient(X1):
    S1=-symToVecCalc(fg,X1); %= -fg(X1)
    S1=S1./max(abs(S1));
    %S1=S1./(sqrt(sum(S1.^2)));
    normalize_s=false;
    %%
    %Get Lambda Star (Optimal Step Size):
    try %Cubic:
        if(Fletcher_Reeves_params.quad_interp_only)
               error('Quad Interpolation Only. Go to catch phrase.');
        end
        [lambda_star,accuracy] = SymbCubicInterpOpt(f,fg,t,X1, S1, epsilon1, epsilon2,lmax,normalize_s,verbose);
    catch %With Quadratic
        [lambda_star,accuracy] = SymbQuadInterpOpt(f,t,X1,S1,epsilon1,normalize_s,verbose);
    end

    X2= X1+lambda_star*(S1); %Remember: S1 was normalized.

    %============= Prepare for the loop: ==================
    prev_X=X1;
    prev_grad=symToVecCalc(fg,prev_X);
    mag_prev_grad_sq=sum(prev_grad.^2);
    prev_f=symToVecCalc(f,prev_X);

    current_X=X2;
    current_grad=symToVecCalc(fg,current_X);
    mag_current_grad_sq=sum(current_grad.^2);
    current_f=symToVecCalc(f,current_X);

    prev_S=S1;

    exit_loop=false;
    %%%Eof Loop preparations ===============================

    % =========================== MAIIN LOOP: =================================
    num_iterations=0;
    while(~exit_loop)
       num_iterations=num_iterations+1; %Count iterations (One iteration was already done before)

       beta_star=(mag_current_grad_sq/mag_prev_grad_sq);
       current_S= -current_grad + beta_star.*prev_S;
       %current_S= current_S./max(abs(current_S)); Note: NO NORMALIZATION EVEN IN CUBIC FUNCTION !

       %Get Lambda Star (Optimal Step Size):
       try %Cubic
           if(Fletcher_Reeves_params.quad_interp_only)
               error('Quad Interpolation Only. Go to catch phrase.');
           end
           [lambda_star,accuracy] = SymbCubicInterpOpt(f,fg,t,current_X,current_S, epsilon1, epsilon2,lmax,normalize_s,verbose);
       catch %With Quadratic
           [lambda_star,accuracy] = SymbQuadInterpOpt(f,t,current_X,current_S,epsilon1,normalize_s,verbose);
       end
        
       
       if(prev_f>current_f) %Change to quadratic.. cubic is no longer working
           Fletcher_Reeves_params.quad_interp_only=true;
       end
       prev_X=current_X;
       prev_grad=current_grad;          
       mag_prev_grad_sq=mag_current_grad_sq;
       prev_f=current_f;


       current_X= prev_X+lambda_star.*current_S;
       current_grad=symToVecCalc(fg,current_X);
       mag_current_grad_sq=sum(current_grad.^2);
       current_f=symToVecCalc(f,current_X);

        prev_S=current_S;
        %=================Optimality Check:===================
        %Debugging : Print the magnitude of the gradient
        %Second condition (On Gradient): 
        if(mod(num_iterations,Fletcher_Reeves_params.print_every)==0)
            disp(current_f)
        end
        
        
        
        %Equivalent to ORing all conditions:
        %First condition (On Function): 
        %if( abs((current_f - prev_f)  / current_f)  <= eps1 )
        %    exit_loop=true;
        %    break;
        %end

        
        if( sqrt(sum(current_grad.^2))  <= eps2 )
            exit_loop=true;
            break;
        end

        %Third condition (On X):
        %if( (sqrt( sum( (current_X - prev_X).^2 ) )  / sqrt( sum(current_X.^2) )) <= eps3 )
        %   exit_loop=true;
        %   break;
        %end

        %Eof Optimality Check ================================


    end

    %==Eof: MAIIN LOOP: =======================================================

    f_star=current_f;
    X_star=current_X;


    
end