function [index] = findFirstOccurNotIn(A,x,forbiddenIndices)
    %Find index of first occurence for value x in vector A, given the index
    %is not in forbiddenIndices
    index=-1; %Returns -1 if value x not found in A
    for i=1:length(A)
        if(A(i)==x)
            if(length(forbiddenIndices(forbiddenIndices==i))==0)
                index=i;
                break;
            end
        end
    end
    


end