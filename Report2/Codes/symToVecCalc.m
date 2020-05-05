function out = symToVecCalc(f,X)
    
    out=double(subs(f,symvar(f),X'));
    %{
    str_to_eval="f(";
    for i = 1:length(X)
        str_to_eval=strcat(str_to_eval,num2str(X(i)));
        if(i<length(X))
            str_to_eval=strcat(str_to_eval,",");
        end
    end
    str_to_eval=strcat(str_to_eval,")")
    out=double(eval(str_to_eval)); 
    %}

end