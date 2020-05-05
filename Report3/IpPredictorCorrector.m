function [x_star,f_star,history] = IpPredictorCorrector(problem)

    %Nocedal and Wright, algorithm 14.3
    %========================== Problem Paramteres ========================
    if(problem.maximize)
        maximize=1;                 % 1 if maximization, -1 if minimization
    else
        maximize=-1;                %i.e: minimize
    end
    c=problem.c; orig_c=c; c= -1*maximize*c;  % Weights of x's in objective function %Reverse sign if necessary
    A=problem.A;                    % m * n matrix. Weights of x's in m constraint equality equations
    b=problem.b;                    % m*1 vector (column vector)
    sigma=problem.sigma;            % [0:1] (Centering parameter)

    %Problem paramters inside the function:
    m=size(A,1); %number of constraints;
    n=size(A,2); %number of variables
    e=ones(n,1);

    %Initialized Parameters:
    try
        x=problem.init_x; %Edit: Replace with good initialization
    catch
        x=ones(n,1);
    end 


    try
        s=problem.init_s; %Edit: Replace with good initialization
    catch
        s=ones(n,1);
    end   


    try
        lambda=problem.init_lambda; %Edit: Replace with good initialization
    catch
        lambda=zeros(m,1); %Edit: Replace with good initialization
    end   

    X=diag(x);
    S=diag(s);



    %== Eof: Problem Paramteres ============================================


    %Check (x,s) >=0:
    assert(all(x>0),"Error.. the input x has <= 0  elements..Exiting")
    assert(all(s>0),"Error.. the s has <= 0 elements.. Exiting")


    %Prepare for loop:
    iters=0;
    mu= ( x' * s )  ./n ;

    %Store History:
    x_hist=[x];
    s_hist=[s];
    mu_hist=[mu];


    while( mu >problem.StopMu)
        %F:
        rb= A*x - b;
        rc= A'*lambda + s - c;
        F=[rc;rb;X*S*e];

        %R_(2n+m) matrix (Jacboian):: J :
        %{  
            --                                --     
            |  0_(n,n)    AT(n,m)     I(n,n)    |
            |                                   |
            |   A(m,n)    0_(m,m)     0(m,n)    |
            |                                   |
            |   S(n,n)    0_(n,m)     X(n,n)    |
            --                                 --

        %}
        J=[
            zeros(n,n)  ,   A'          ,  eye(n)       ; ...
            A           ,  zeros(m,m)   ,  zeros(m,n)   ; ...
            S           ,  zeros(n,m)   ,  X            ; ...
            ];

        %   Deltas Vector for Affine transformations :
        D_aff=J\(-F); %Inefficient, but works.
        dx_aff      =    D_aff(1:n,:);
        dlambda_aff =    D_aff(n+1:n+m,:);
        ds_aff      =    D_aff(n+1+m:end,:);


        %   Take step:
        eps=1e-12;
        dxdes=dx_aff(dx_aff<0);
        xdes=x(dx_aff<0);

        dsdes=ds_aff(ds_aff<0);
        sdes=s(ds_aff<0);


        alpha_aff_pri = min([1;(eps+problem.gamma*mu-xdes)./dxdes]);
        alpha_aff_dual= min([1;(eps+problem.gamma*mu-sdes)./dsdes]);



        %Don't actually move a step, just get mu_aff:
        mu_aff = ( x + alpha_aff_pri.* dx_aff)' * ( s + alpha_aff_dual.* ds_aff ) /n;

        sigma  = (mu_aff/mu)^3;

        if(sigma>=1)
            sigma=0.99; %TO avoid infinite loop
        elseif(sigma<0)
            sigma=0001;
        end

        dXaff=diag(dx_aff);
        dSaff=diag(ds_aff);


        %Calculate corrected Deltas vector:
        D=J\( -F + [zeros(size(rc));zeros(size(rb));(mu*sigma.*e - dXaff*dSaff*e)] );


        %Take Actual Step:
        %==================================================================%
        dx      =    D(1:n,:);
        dlambda =    D(n+1:n+m,:);
        ds      =    D(n+1+m:end,:);


        eps=1e-12;
        dxdes=dx(dx<0);
        xdes=x(dx<0);

        dsdes=ds(ds<0);
        sdes=s(ds<0);




        alpha_true_pri = min([1;(eps-xdes)./dxdes]);
        alpha_true_dual= min([1;(eps-sdes)./dsdes]);

        x       =      x       +   alpha_true_pri.* dx       ;
        s       =      s       +   alpha_true_dual.* ds      ;
        lambda  =    lambda    +   alpha_true_dual.*dlambda  ;

        X=diag(x);
        S=diag(s);
        mu= (x'*s) /n;
        
        %Store history:
        x_hist=[x_hist,x];
        s_hist=[s_hist,s];
        mu_hist=[mu_hist,mu];
        %==================================================================%


        %Check (x,s) >0
        assert(all(x>0),"Error.. the input x has <= 0  elements..Exiting")
        assert(all(s>0),"Error.. the s has <= 0 elements.. Exiting")

        iters=iters+1;
    end

    iters=iters+1; %To account for initial state

    history=struct();
    history.x_hist=x_hist;
    history.s_hist=s_hist;
    history.mu_hist=mu_hist;
    history.iters=iters;

    x_star=x;
    f_star=orig_c'*x;
    
end