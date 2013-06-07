function AAs = MoveEqualVar(AAs,orig,copy,format)

switch format
    case 'transp'
        for i=1:length(AAs)
            A = AAs{i}';
            tmp = A(:,copy(1));
            for j=2:length(copy)
                tmp = tmp + A(:,copy(j));
            end
            A(:,orig) =A(:,orig) + tmp;
            AAs{i} = A';
            
        end
    case 'normal'
        for i=1:length(AAs)
            tmp = AAs{i}(:,copy(1));
            for j=2:length(copy)
                tmp = tmp + AAs{i}(:,copy(j));
            end
            AAs{i}(:,orig) = AAs{i}(:,orig) + tmp;
        end
        
end
