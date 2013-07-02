function full_name = BLOM_GetPolyBlockEquation(A,C,gA,gC,names,max_rows,max_length)

full_name = [];
for i = 1:min(size(C,1),max_rows)
    fname = CreatePolyFromMatrix(A,C(i,:),names,'low');
    if (length(fname) > max_length-5)
        fname = 'f(x)';
    end
    
    if (gC(1) == 1 && size(gA,1)==1 && sum(abs(gA(1,1:end)))==0)
        name = sprintf('%s',fname);
    else
        if size(gC,1) == 1
            gname = CreatePolyFromMatrix(gA,gC,names,'low');
        else
            gname = CreatePolyFromMatrix(gA,gC(i,:),names,'low');
        end
        name = sprintf('(%s)/(%s)',fname,gname);
        if (length(name) > max_length)
            
            gname = 'g(x)';
            name = sprintf('(%s)/%s',fname,gname);
        end
    end
    full_name = [full_name   name  char(10)];
end

