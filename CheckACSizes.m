function CheckACSizes(A,C,gA,gC)

if (size(A,2) ~= size(gA,2))
    str =[gcb ': Number of columns of A in f must be equal to number of columns of A in g!']; 
    disp(str);
    error(str);
end

if (size(A,1) ~= size(C,2))
    str =[gcb ': Number of rows of A in f must be equal to number of columns of C in f!'];
    disp(str);
    error(str);
end

if (size(gA,1) ~= size(gC,2))
    str =[gcb ': Number of rows of A in g must be equal to number of columns of C in g!'];
    disp(str);
    error(str);
end

if (size(gC,1) ~= 1)
    str =[gcb ': g must be a scalar function'];
    disp(str);
    error(str);
end
