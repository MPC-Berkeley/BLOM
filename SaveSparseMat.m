function SaveSparseMat(data,fname)
[i,j,val] = find(data);
if (size(i,1) == 1)
    data_dump = [i;j;val];
else
    data_dump = [i,j,val]';
end
if isempty(data) || (data(end,end) == 0) % no last element
    data_dump = [data_dump [size(data,1); size(data,2); 0]];
end

f = fopen(fname,'wt');
fprintf(f,'%d %d %.18g\n',data_dump);
fclose(f);