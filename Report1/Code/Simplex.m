function [x_star, f_star, num_iterations, history] = Simplex(lin_problem)



    %======================I N P U T S==========================%

    c = lin_problem.c;                          % 1*n vector (NOTE: note like course slides where c = n * 1)
    A = lin_problem.A;                          % Coefficients of x in inequalities.. size= m * n
    b = lin_problem.b;                          % Right hand side of inequalitise

    x=[1:length(c)]';                           % Note the tanspose (column vector)
    xs=[x(end)+1:x(end)+size(A,1)]';              % Note the tanspose (column vector)
    xvals=zeros((length(x)+length(xs)),1);      % Initialized with zeros.. to store xvalues

    %If maximimze, take -c:
    if(lin_problem.maximize)
        c=-1 *[c,zeros(1,length(xs))];
    else
         c=[c,zeros(1,length(xs))];
    end

    %Aggregate A:
    A=[A,eye(length(xs))];              
    %Aggregate (useful  but redundant data):
    xall_inds=[x;xs];
    %==== Eof: I N P U T S=====================================%



    %Initially:
    xB=xs; %Init xB 
    B=A(:,xs'); %init B %remember to Update back B into A
    Binv=inv(B); %Initially = B

    num_iterations=0;
    
    history   = struct();
    history.z = [];
    while(true)
        %Step:
        xvals(xB)= Binv*b;
        cB= c(xB);
        z= cB*xvals(xB);
        history.z=[history.z,z];


        %2. Determine entering basic variable:
        %Entering= most negative coefficient
        %Get Value then index of most negative non-basic co-eff
        c_entering_val=min(c(removeSharedEl(xall_inds,xB)));
        if(c_entering_val>=0)
            break;
            %Edit: Stop here.. optimal
        end
        c_entering_ind=findFirstOccurNotIn(c,c_entering_val,xB);

        %3. Determine leaving non-basic variable
        %Get co-effs of entering variable in each equation from 
        %Get coefficients of entering variable from B:
        coeffs=A(:,c_entering_ind);
        %Get indices of strictly positive coefficients (if not, stop, infeasible):
        positive_coeef_inds=find(coeffs>0);
        if(length(positive_coeef_inds)==0)
           %Infeasible.. stop and exit 
        end
        %      Find minimum non-negative ratio value and index:

        temp_rhs=xvals(xB);
        pos_coeffs_ratios = temp_rhs(positive_coeef_inds)./ coeffs (positive_coeef_inds);
        pos_coeffs_ratios(pos_coeffs_ratios<0)=max(pos_coeffs_ratios)+1; %Change negative values to large positive values so that the min function ignores them
                                                                         %Has no use in the algorithm semantics but eliminates negative ratios from consideration
        [min_pos_ratio_value,ind_in_positive_coeef_inds]=min(pos_coeffs_ratios);

        index_of_min_non_negative_ratio_in_A_rows=positive_coeef_inds(ind_in_positive_coeef_inds);
        leaving_index=index_of_min_non_negative_ratio_in_A_rows;

        %Now we got the pivot column in (leaving_index,c_entering_ind) --> (row, column)

        xB(leaving_index)= c_entering_ind;

        B=A(:,xB); %init B %remember to Update back B into A
        Binv=inv(B); %Initially = B
        
        num_iterations=num_iterations+1;
    end





    x_star=xvals(x);
    if(lin_problem.maximize)
        history.z= - history.z;
        f_star= -1 * z;
    else
        f_star = z;
    end
    
    %num_iterations;
    





end