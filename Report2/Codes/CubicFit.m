%Fit Cubic function given [A,B], f(A),f'(A),f(B),f'(B). 

function dcba = CubicFit(x,y,verbose)
    %Returns a cubic equation fit for exactly 2 points and the deirivative of the function at those (
    % two points (Analytically). 
    %The input points are named A,B (capital). The quadratic coefficients are (d,c,b,a)
    %for dx^3 +cx^2+bx+a #NOTE THE NAMING CONVENTION !.
    %returned [d,c,b,a]
    %%%=========================INPUTS================================%
    %x: array of two inputs points
    %y: array of four elements = f(A),f'(A),f(B),f'(B). Where ' denotes
        %derivative.
    %verbose: true to display output, false to not display output. 
    %%%====Eof: INPUTS=================================================%

    assert(length(x)==2,"Error: Please note that the length of x should be 2.\n In the form: [A,B] \nExiting..");
    assert(length(y)==4,"Error: Please note that the length of y should be 4.\n In the form: [f(A), f'(A), f(B), f'(B)]. \nExiting..");

    A=x(1);
    B=x(2);
    fA=y(1);
    fdA=y(2);
    fB=y(3);
    fdB=y(4);

    assert(A~=B, "Error: Your two selected points must not be equal.");

    Z=( 3*(fA-fB) / B )+ fdA + fdB; %BLock
    d= (1 / (3*(A-B)^2) ) * ( 2*Z + fdA + fdB );
    c= -(1/(A-B)^2) * ...
        ( (A+B)*Z + B*fdA + A*fdB );
    b= (1/(A-B)^2) * ...
       ( B^2*fdA + A^2*fdB + 2*A*B*Z );
    a= fA - b*A - c*A^2 - d*A^3;

    if(verbose)
        fprintf("%.2fx^3 + %.2fx^2 + %.2fx + %.2f \n",d,c,b,a); 
    end

    dcba=[d,c,b,a];
end