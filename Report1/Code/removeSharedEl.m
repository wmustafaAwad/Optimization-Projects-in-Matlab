function [A] = removeSharedEl(A,B)
    %Removes elements from B found in A
    for i=1:length(B)
        A(A==B(i))=[];
    end
end