function [x_star,f_star,history] = IpPrimDualPathFollowing(problem)
    %Nocedal and Wright, algorithm 14.1
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
    x=ones(n,1); %Edit: Replace with good initialization
    X=diag(x);

    s=ones(n,1); %Edit: Replace with good initialization
    S=diag(s);

    lambda=zeros(m,1); %Edit: Replace with good initialization

    %== Eof: Problem Paramteres ============================================


    %Check (x,s) >=0:
    assert(all(x>0),'Error.. the input x has <= 0  elements.. Exiting')
    assert(all(s>0),'Error.. the s has <= 0  elements.. Exiting')


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

        %   Deltas Matrix:
        D=J\(-F + [zeros(size(rc));zeros(size(rb));(mu*sigma.*e)]); %Inefficient, but works.
        dx=D(1:n,:);
        dlambda=D(n+1:n+m,:);
        ds= D(n+1+m:end,:);


        %   Take step:
        eps=1e-13;

        dxdes=dx(dx<0);
        xdes=x(dx<0);

        dsdes=ds(ds<0);
        sdes=s(ds<0);


        alpha_pri = min([1;(eps-xdes)./dxdes]);
        alpha_dual= min([1;(eps-sdes)./dsdes]);
        x       =      x       +   alpha_pri.* dx      ;
        s       =      s       +   alpha_dual.* ds      ;
        lambda  =    lambda    +   max([alpha_pri,alpha_dual]).*dlambda  ;
        
        X=diag(x);
        S=diag(s);
        mu= (x'*s )/n;

        %Store history:
        x_hist=[x_hist,x];
        s_hist=[s_hist,s];
        mu_hist=[mu_hist,mu];



        %Check (x,s) >0
        assert(all(x>0),'Something went wrong.. the x has <= 0  elements')
        assert(all(s>0),'Something went wrong.. the s has <= 0  elements')

        iters=iters+1;
    end
    
    iters=iters+1;


    history=struct();
    history.x_hist=x_hist;
    history.s_hist=s_hist;
    history.mu_hist=mu_hist;
    history.iters=iters;

    x_star=x;
    f_star=orig_c'*x;
    
end