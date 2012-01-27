function CreateIpoptDAT(name,fixed,x0)

f = fopen([ name 'X0.dat'],'wt');
fprintf(f,'%16.15e \n',x0);
fclose(f);

f = fopen([ name 'Fixed.dat'],'wt');

Data = zeros(length(fixed.AAs),1);
for i=1:length(fixed.AAs)
    Data(i) = fixed.Cs{i}(2);
end

fprintf(f,'%16.15e \n',Data);
fclose(f);
