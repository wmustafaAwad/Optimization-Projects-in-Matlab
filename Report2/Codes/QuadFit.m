%Fit quadratic function given three points A,B,C & f(A),f(B),f(C):
function cba = QuadFit(x,y,verbose)
%Returns a quadratic equation fit for exactly 3 points  (Analytically). The
%points are named A,B,C (capital). The quadratic coefficients are (c,b,a)
%for cx^2+bx+a #####: NOTE THE NAMING CONVENTION !
%returned [c,b,a]
%%%=========================INPUTS================================%
%x: array of three inputs points
%y: array of three elements = f(x)
%verbose: true to display output, false to not display output. 
%%%====Eof: INPUTS=================================================%


%{
%%Sample Run:
verbose=true;
x=[0,1,2];
y=x.^2;
QuadFit(x,y,verbose)

=============================================================
Output:
1.00x^2 + 0.00x + 0.00 

ans =

     1     0     0


%}

assert(length(x)==3 && length(y)==3,"Error: Please note that either the length of x, the length of y, or both; is not 3. Exiting..");

A=x(1);
B=x(2);
C=x(3);
fA=y(1);
fB=y(2);
fC=y(3);

assert(A~=B && B~=C && A~=C, "Error: All of your three selected points must not be equal. At least two of them are equal.");

%Analytically Derived Equations:
b= ( fA*(B^2-C^2)  +  fB*(C^2-A^2) + fC*(A^2-B^2) )/ ...
    ( A*(B^2-C^2)  +  B*(C^2- A^2) + C* (A^2-B^2) );
c=( fB-fC- b*(B-C) )/ (B^2-C^2);
a=( fA+fB -b*(A+B) - c*(A^2+B^2) )/2;

if(verbose)%Print if asked to
    fprintf("%.2fx^2 + %.2fx + %.2f \n",c,b,a); 
end
cba=[c,b,a];
end