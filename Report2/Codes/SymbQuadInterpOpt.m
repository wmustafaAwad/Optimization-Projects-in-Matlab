function [lambda_star,accuracy,n_refits] = SymbQuadInterpOpt(f,t,x,s,epsilon1,normalize_s,verbose)
    %Symbolic Quadratic Interpolation for 1-D Optimization
    %Given INPUTS, find the best lambda to minimize the  function f,x: input,s:direction, epsilon1: condition
    %on on function approximation.
    %Source: Rao, Engineering Optimization. Non-linear Programming chapter.

    %%%=========================INPUTS================================%
    %f; Symbolic function of x1,x2...xn (symbolic variables are defined
    %outisde SymbQuadInterpOpt
    %t; Starting step; Initially: B-A=t, C-B=t; C-A=2*t;
    %epsilon1; %Required Accruacy. judged by:
        %abs((f_lambda_star-h_lambda_star)/f_lambda_star)<=epsilon1
    %x; Input (initial) point (has dimensions= expected input by funct_name's function
    %s; Direction in which to search for optimum lambda. length(s) must be equal
        %to length(x).
    %verbose; %If true: print all quadratic fits. If false: print nothing.
    %%%====Eof: INPUTS=================================================%
    if(normalize_s)
        s=s./max(s); %Normalize direction.
    end
    %========== Step 1 ==========%
    A=0;
    fA=symToVecCalc(f,x+A.*s);%funct(x+A.*s);
    f1=symToVecCalc(f,x+(A+t).*s);
    exit_loop=false;
    n_refits=0;
    %== Eof: Step 1 =============%


    while(~exit_loop)
        %First Loop (Steps 2-5)
        if(f1>fA)
            fC=f1;
            C=A+t; %A is initially zero, but not necessarilty zero later
            B=A+t/2;
            fB=symToVecCalc(f,x+B.*s);
            t=t/2;%%ADDED, nothing changed
            exit_loop=true;
        else%if(f1<=fA)
            fB=f1;
            B=A+t;
            f2=symToVecCalc(f,x+(A+2*t).*s);%A is initially zero, but not necessarilty zero later
            if(f2>=f1)
                fC=f2;
                C=A+2*t;
                exit_loop=true;
            else%if(f2<f1)
                f1=f2;
                t=2*t;
            end
        end
    end
    %==============================================%
    %Calculate lambda_star: (Part of steps 2-5)
    cba=QuadFit([A,B,C],[fA,fB,fC],verbose);
    c=cba(1);
    b= cba(2);
    a=cba(3);
    lambda_star= -b/(2*c); %A-b/(2*c)
    h_lambda_star=c*lambda_star^2+b*lambda_star+a;
    f_lambda_star=symToVecCalc(f,x+lambda_star.*s);
    %==============================================%
    %Check if we are done: (Do we need to refit ?)
    if(abs((f_lambda_star-h_lambda_star)/f_lambda_star)<=epsilon1)
        %Keep exit_loop=true
    else
        exit_loop=false; %Now we need to go into refitting loop
    end
    %==============================================%
    while(~exit_loop)
       %So, we need to refit (second loop--Table 5.5):
       n_refits=n_refits+1;
       %TRY REFITTING: (Four cases)
       if(lambda_star>=B)
           if(f_lambda_star<=fB)
               Anew=B;
               fAnew=fB;
               %
               Bnew=lambda_star;
               fBnew=f_lambda_star;
               %
               Cnew=C;
               fCnew=fC;
           else%f_lambda_star>fB
               Anew=A;
               fAnew=fA;
               %
               Bnew=B;
               fBnew=fB;
               %
               Cnew=lambda_star;
               fCnew=f_lambda_star;
           end
       else%lambda_star<B
           if(f_lambda_star<=fB)
               Anew=A;
               fAnew=fA;
               %
               Bnew=lambda_star;
               fBnew=f_lambda_star;
               %
               Cnew=B;
               fCnew=fB;
           else%f_lambda_star>fB
               Anew=lambda_star;
               fAnew=f_lambda_star;
               %
               Bnew=B;
               fBnew=fB;
               %
               Cnew=C;
               fCnew=fC;
           end
       end
       %Now, test if refitting has worked:
       cba= QuadFit([Anew,Bnew,Cnew],[fAnew,fBnew,fCnew],verbose);
       c=cba(1);
       b= cba(2);
       a=cba(3);
       lambda_star= -b/(2*c);
       h_lambda_star=c*lambda_star^2+b*lambda_star+a;
       f_lambda_star=symToVecCalc(f,x+lambda_star.*s);
       if(abs((f_lambda_star-h_lambda_star)/f_lambda_star)<=epsilon1)
           %then we have got a good lambda star, leave
           exit_loop=true;
       else %Still, after refitting ,we have got a bad lambda
           %Just update, and don't exit (remember, exit_loop is false here). Go back to step 2
           A=Anew;
           fA=fAnew;
           %
           B=Bnew;
           fB=fBnew;
           %
           C=Cnew;
           fC=fCnew;
       end

    end

    accuracy=abs((f_lambda_star-h_lambda_star)/f_lambda_star);

end