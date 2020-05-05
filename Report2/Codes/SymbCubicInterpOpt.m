function [lambda_star,accuracy,n_refits] = SymbCubicInterpOpt(f,fg,t,X, S, epsilon1, epsilon2,lmax,normalize_s,verbose)
    %Version of CubicInterpOpt that has input function as symbolic 
    %============================INPUTS=================================%
    %f: function of x1,x2,x3...xn (all are symbolic) and defined outside
    %the function
    %fg: gradient of f (symbolic column-vector)
    
    %t: starting step. B-A=t
    %X: point to start from (length(X) = n)
    %S: search direction. The normalization method inside the function for
    %S is S= S=S./max(abs(S));. Remember to use the normalized S with the direction

    %epsilon1: accuracy constraint on function
    %epsilon2: (not currently used) accuracy constraint on function derivative
    %lmax: maximum number of refits
    
    %verbose: if true, print all cubic fitting function coefficients

    %%%====Eof: INPUTS=================================================%
    if(normalize_s)
        S=S./max(abs(S));
    end
    %%%====Eof: INPUTS=================================================%
    %%
    %A:
    A=0;
    fA=symToVecCalc(f,X+A.*S); %f(lambda)
    fdA=S'*symToVecCalc(fg,X+A.*S); %df/d(lambda)
    if(fdA>0)
        error('Chosen direction is not a descent direction. Exiting..\n')
    end
    %%
    %B:
    exit_loop=false;
    loop_cntr=0;
    while(~exit_loop)
       fB=symToVecCalc(f,X+t.*S); %double(f(t));
       fdB=S'*symToVecCalc(fg,X+t.*S); %double(fd(t));
       if(fdB>0)
           l=0;
           B=t;
           exit_loop=true;
       else
           %SpeedUp Convergence:
           %fA=fB;
           %fdA=fdB;
           %A=t;
           %----%
           t=2*t;
       end
       loop_cntr=loop_cntr+1;
       if(loop_cntr>100 || abs(fdB)>1e20)
           error("Cubic Interpolation failed to retrieve a second point B. Please check your f(lambda) is not quadratic, and that direction input is a descent direction.");
       end
    end

    exit_loop=false; %New loop to start
    n_refits=0;
    while(~exit_loop)
        n_refits=n_refits+1;
        Z= (3*(fA-fB)/(B-A)) + fdA + fdB ;
        Q= sqrt(Z^2 - fdA*fdB);
        if(~isreal(Q))
            warning("Q is not real for some reason. This means that f'(A) and f'(B) do not have opposite signs.\n Please make sure that your input direction is a descent direction.\n")
        end
        lambda_star_p= A+ (fdA + Z + Q)*(B-A)/(fdA+fdB+ 2* Z); %plus Q
        lambda_star_m= A+ (fdA + Z - Q)*(B-A)/(fdA+fdB+ 2* Z); %minus Q
        f_lambda_star_p= symToVecCalc(f,X+lambda_star_p.*S); %double(f(lambda_star_p));
        f_lambda_star_m= symToVecCalc(f,X+lambda_star_m.*S); %double(f(lambda_star_m));


        if(f_lambda_star_p<f_lambda_star_m)
            lambda_star=lambda_star_p;
            f_lambda_star=f_lambda_star_p;
        else %use the other lambda
            lambda_star=lambda_star_m;
            f_lambda_star=f_lambda_star_m;
        end

        

        l=l+1;

        %Fit:
        dcba= CubicFit([A,B],[fA,fdA,fB,fdB],verbose);
        d= dcba(1);
        c= dcba(2);
        b= dcba(3);
        a= dcba(4);
        h_lambda_star=d*lambda_star^3 + c*lambda_star^2 ...
            + b*lambda_star + a;

        if(abs((h_lambda_star-f_lambda_star)/f_lambda_star)<=epsilon1 && S'*symToVecCalc(fg,X+lambda_star.*S)<=epsilon2)
            exit_loop=true; %take current lambda_star

        else %First convergence condition not satisfied
            if(l>=lmax) %stop, that's enough
                error('Max Number of Allowed iterations Exceeded');
                exit_loop=true; %take current lambda_star
            else %choose between A,B,lambda_star
                fd_lambda_star=S'*symToVecCalc(fg,X+lambda_star.*S);%double(fd(lambda_star));
                if(fd_lambda_star<0)%::> increasing lambda decreases function
                    A=lambda_star;
                    fA=f_lambda_star;
                    fdA=fd_lambda_star;
                    %B is left as is
                else %fd_lambda_star>=0::> increasing lambda increases function
                    B=lambda_star;
                    fB=f_lambda_star;
                    fdB=fd_lambda_star;
                    %A is left as is
                end
            end
        end


    end

    accuracy=abs((h_lambda_star-f_lambda_star)/f_lambda_star);
end