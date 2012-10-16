function str = CreatePolyFromMatrix(A,C,names,precision,target)

if (nargin < 4)
    precision = 'high';
end

if (nargin < 5)
    target = 'matlab';
end

str = '';

idx = find(C);
if (~isempty(idx))
    for i = idx % 1:size(A,1)
        %     if (C(i) == 0)
        %         continue;
        %     end
        [jA ] = find(A(i,:));
        skip_mult = false;
        if C(i)==1 &&  ~isempty(jA)
            skip_mult = true;
            if (i>1)
                str = [str  '+'];
            end
        elseif C(i)==-1 &&  ~isempty(jA)
            skip_mult = true;
            str = [str  '-'];
        elseif (C(i)>0 && (i > 1))
            switch precision
                case 'high'
                    str = [str  '+'   sprintf('%15.14g',full(C(i)))];
                otherwise
                    str = [str  '+'   num2str(full(C(i)))];
            end
            
        else
            switch precision
                case 'high'
                    str = [str     sprintf('%15.14g',full(C(i)))];
                otherwise
                    str = [str     num2str(full(C(i)))];
            end
        end
        
        
        
        
        for j = 1:length(jA)
            if (A(i,jA(j)) > 0 && mod(A(i,jA(j)),1) == 0 && A(i,jA(j)) <  min(BLOM_FunctionCode)    ) % integer power
                for k = 1:A(i,jA(j))
                    if (~(skip_mult && j ==1 && k ==1) )
                        str = [str  '*'];
                    end
                    str = [str  names{jA(j)} ];
                end
            else      % non integer power, negative or exp/log
                if (~(skip_mult && j ==1 ) )
                    str = [str  '*'];
                end
                if any(A(i,jA(j)) == BLOM_FunctionCode)
                    if A(i,jA(j)) == BLOM_FunctionCode('exp')
                        str = [str  'exp(' names{jA(j)} ')' ];
                    elseif A(i,jA(j)) == BLOM_FunctionCode('log')
                        str = [str  'log(' names{jA(j)} ')' ];
                    elseif A(i,jA(j)) == BLOM_FunctionCode('sin')
                        str = [str  'sin(' names{jA(j)} ')' ];
                    elseif A(i,jA(j)) == BLOM_FunctionCode('cos')
                        str = [str  'cos(' names{jA(j)} ')' ];
                    elseif A(i,jA(j)) == BLOM_FunctionCode('tanh')
                        str = [str  'tanh(' names{jA(j)} ')' ];
                    else
                        error(['Unrecognized function code ' num2str(A(i,jA(j)))])
                    end
                else % non integer power or negative
                    switch(target)
                        case 'matlab'
                            switch precision
                                case 'high'
                                    str = [str names{jA(j)} '^' sprintf('%15.14g',full(A(i,jA(j))))];
                                otherwise
                                    str = [str names{jA(j)} '^' num2str(full(A(i,jA(j))))];
                            end
                        case 'C'
                            str = [str 'pow('  names{jA(j)} ',' sprintf('%15.14g',full(A(i,jA(j)))) ')' ];
                            
                    end
                end
                
            end
            
        end
        
    end
end
if isempty(str)
    str = '0';
end