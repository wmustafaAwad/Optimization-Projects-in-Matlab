function [X_star,f_star,i] = GoldenOpt(f,x,s,golden_params)
%returns min_pair: (x,f(x)) that are either minimum or maximum of the input function using the 0.618 (golden) ratio.
%%%=========================INPUTS================================%
    %xl=region(1);               %lower X in search region
    %xu=region(2);               %uppwer X in search region
    %accuracy;                   %region in which to search for optimal pair
    %minimize;                   %If true: find minimizing pair, if false: find maximizing pair
    %temp_f=str2func(funct_name);%funct_name: name of the function to be used for minimzation
    %verbose;                    %if true: print each step, if false: just return output
%%%====Eof: INPUTS=================================================%
    %Golden Section
    %{
    %%Sample Run:
    region=[0,4];
    accuracy=1e-2;
    minimize=false;
    funct_name="temp_f";
    verbose=false;
    GoldenOpt(region,accuracy,minimize,funct_name,verbose)
    returns: 1.442719099991563   1.775474795224482 
    %%Those are: max_x, max_f(x)
    %}
    
    %1- Given xl,xu initial, r=1/1.61803398875, accuracy, num_steps=(0.618)^n such that < desired accuracy.
    
    %%%=========================INPUTS================================%
    xl=golden_params.region(1); %lower X in search region
    xu=golden_params.region(2); %uppwer X in search region
    accuracy=golden_params.eps;    %region in which to search for optimal pair
    minimize=golden_params.minimize;    %If true: find minimizing pair, if false: find maximizing pair
    temp_f=f; %funct_name: name of the function to be used for minimzation
    verbose=golden_params.verbose; %if true: print each step, if false: just return output
    %%%====Eof: INPUTS=================================================%
    
    %%%=======================Initializations===========================%
    f_arr=[nan,nan,nan,nan]; %f_xl, f_x2, f_x1,f_xu
    f_arr_new=[nan,nan,nan,nan];
    x_arr=[nan,nan,nan,nan];%xl, x2, x1, xu
    x_arr_new=[nan,nan,nan,nan];
    r=1/1.61803398875;
    num_steps=ceil(log(accuracy)/log(r)); %%For Golden Section Only
    %2- d=r*(xu-xl)
    d=r*(xu-xl);
    d_new=d;
    x_arr(1)=xl;
    x_arr(4)=xu;
    %3- x1=xl+d
    x_arr(3)=xl+d; %x_arr(3)=x_1
    %4- x2=xu-d
    x_arr(2)=xu-d; %x_arr(2)=x2
    %5- Calculate (or retrieve) %f_xl, f_x2, f_x1,f_xu
    f_arr(1)=symToVecCalc(temp_f,x+x_arr(1).*s);%temp_f(x_arr(1));
    f_arr(2)=symToVecCalc(temp_f,x+x_arr(2).*s);%temp_f(x_arr(2));
    f_arr(3)=symToVecCalc(temp_f,x+x_arr(1).*s);%temp_f(x_arr(3));
    f_arr(4)=symToVecCalc(temp_f,x+x_arr(4).*s);%temp_f(x_arr(4));
    X_star=x_arr(2);
    f_star=f_arr(2);
    %%%====Eof: Initializations========================================%
    
    
    
    
    
    
    %FOR PRINTING ONLY:
    f_arr_new=f_arr;
    x_arr_new=x_arr;
    if(verbose)
        fprintf(['\n---------------------------------------------------------------------------------------------------\n', ...
            'i:%d|xl:%.4f|f(xl):%.4f|x2:%.4f|f(x2):%.4f|x1:%.4f|f(x1):%.4f|xu:%.4f|f(xu):%.4f|d:%.4f\n', ...
            '---------------------------------------------------------------------------------------------------'], ...
            1,x_arr_new(1),f_arr_new(1),x_arr_new(2),f_arr_new(2),x_arr_new(3),f_arr_new(3),x_arr_new(4),f_arr_new(4),d_new)
    end
    
   
    %%%================================ L O O P   B E G I N S =================================%%%
    
    for i=2:num_steps
        if((f_arr(3)>=f_arr(2) && minimize) || (f_arr(3)<f_arr(2) && ~minimize)) %if(f(x1)>=f(x2)), take left side %x_arr: %xl, x2, x1, xu
            d_new=r*(x_arr(3)-x_arr(1)); %d'=r*(xu'-xl') = r*(x1-xl) 
            x_arr_new(1)=x_arr(1); %xl'=xl %nothing, remove
            f_arr_new(1)=f_arr(1);
            
            x_arr_new(2)=x_arr(3)-d_new; %x2'=xu'-d'= x1-d'
            f_arr_new(2)=symToVecCalc(temp_f,x+x_arr_new(2).*s);%temp_f(x_arr_new(2)); %ONLY ONE ACTUALLY CALCULATED
            
            
            x_arr_new(3)=x_arr(2); %x1'=x2
            f_arr_new(3)=f_arr(2);
            
            x_arr_new(4)=x_arr(3); %xu'=x1;
            f_arr_new(4)=f_arr(3);
            
        else % f(x1)<f(x2) or not defined %NOTE: EDIT . WILL GO INTO INFINITE LOOP if undefined function! Handle in future versions.
            
            d_new=r*(x_arr(4)-x_arr(2)); %d'=r*(xu'- xl')=r*(xu-x2)
            
            x_arr_new(1)=x_arr(2); %xl'=x2
            f_arr_new(1)=f_arr(2);
            
            x_arr_new(2)=x_arr(3); %x2'=x1,
            f_arr_new(2)=f_arr(3);
            
            x_arr_new(3)=x_arr(2)+d_new; %x1'=xl'+d' = x2+d'
            f_arr_new(3)=symToVecCalc(temp_f,x+x_arr_new(3).*s);%temp_f(x_arr_new(3)); %ONLY ONE ACTUALL CALCULATED
            
            x_arr_new(4)=x_arr(4); %xu'=xu %nothing, remove
            f_arr_new(4)=f_arr(4);
            
        end
        
        
        %Print:
        if(verbose)
            fprintf(['\n---------------------------------------------------------------------------------------------------\n', ...
                'i:%d|xl:%.4f|f(xl):%.4f|x2:%.4f|f(x2):%.4f|x1:%.4f|f(x1):%.4f|xu:%.4f|f(xu):%.4f|d:%.4f\n', ...
                '---------------------------------------------------------------------------------------------------'], ...
                i,x_arr_new(1),f_arr_new(1),x_arr_new(2),f_arr_new(2),x_arr_new(3),f_arr_new(3),x_arr_new(4),f_arr_new(4),d_new)
        end
        
        %Get ready for next loop:
        f_arr=f_arr_new;
        x_arr=x_arr_new;
    end
    %%%==Eof: L O O P ==================================================================================%%%
    
    if(minimize)
        the_min=min(f_arr_new);
        the_min=the_min(1);
        index=find(f_arr_new ==the_min, 1, 'first');
        X_star=x_arr_new(index);
        f_star=f_arr_new(index);
        %min_pair=[x_arr_new(index),f_arr_new(index)];
    else %maximize
        the_min=max(f_arr_new);
        the_min=the_min(1);
        index=find(f_arr_new ==the_min, 1, 'first');
        X_star=x_arr_new(index);
        f_star=f_arr_new(index);
        %min_pair=[x_arr_new(index),f_arr_new(index)];
    end



%{
Steps:
1- Given xl,xu initial, r=/1.61803398875, num_steps=(0.618)^n such that
0.618^n < desired accuracy.
2- d=r*(xu-xl)
3- x1=xl+d
4- x2=xu-d
5- Calculate (or retrieve) f(x1),f(x2),f(xl),f(xu)
Loop for i in range num_steps:
    6- if(f(x1)>=f(x2) and minimize || f(x1)<f(x2) and ~minimize) :: xl'=xl, xu'=x1, x1'=x2 (therefore, f_xl'=f_xl, f_xu'=f_x1, f_x1'=f_x2)
        7- d'=r*(xu'-xl') = r*(x1-xl) 
        8- x2'=xu'-d' , f_x2'=f(x2')   
       else %% f(x1)<f(x2) or not defined:: xl'=x2, xu'= xu, x2'=x1, (therefore, f_xl'=f_x2 f_xu'=f_xu, f_x2'=f_x1)
        7- d'=r*(xu'- xl')
        8- x1'=xl'+d'
    9- Print table (i,xl,f_xl,x2,f_x2,x1,f_x1,xu,f_xu,d)

%}





    
end