function outmatrix=arx_prep(inarray,order)

outmatrix=zeros(order,length(inarray)-order+1);

for i=1:size(outmatrix,2)
    outmatrix(:,i)=inarray(i+order-1:-1:i);
end