function res = CombineFandG(A,C,gA, gC,type)

if (size(gC,1)==1)
    switch(type)
        case 'A'
            res = [A zeros(size(A,1),size(C,1)) ; repmat(gA,size(C,1),1) kron(eye( size(C,1)),ones(size(gA,1),1))  ];
        case 'C'
            res =   [C kron(eye(size(C,1)),-gC)];
    end
elseif size(gC,1)==size(C,1)
    switch(type)
        case 'A'
            res = [A zeros(size(A,1),size(C,1)) ; repmat(gA,size(C,1),1) kron(eye( size(C,1)),ones(size(gA,1),1))  ];
        case 'C'
            res =   [C -gC];
    end
else
    error('number of lines in C matrices of f and g must be either the same, or C of g be a row matrix');
end
% input_type, output_type