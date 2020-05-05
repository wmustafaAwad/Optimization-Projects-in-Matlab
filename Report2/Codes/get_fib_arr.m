function fib_arr = get_fib_arr(n)

% PARAMETERS %
%Returns array of fibonacci numbers from 1 to n (input)
%===================================================================%
%Eg: get_fib_arr(4)
    %returns: [1,1,2,3]
%===================================================================%
    if(n==1 || n==2 || n<1)
        if(n>1)
            fib_arr=[1,1];
        else
            fib_arr=[1];
        end
    else
        fib_arr=zeros(size(1,n));
        fib_arr(1)=1;
        fib_arr(2)=1;
        for i = 3:n
           fib_arr(i)=fib_arr(i-1)+fib_arr(i-2); 
        end
    end
end