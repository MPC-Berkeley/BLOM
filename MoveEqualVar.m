function AAs = MoveEqualVar(AAs,orig,copy)

for i=1:length(AAs)
    tmp = AAs{i}(:,copy(1));
    for j=2:length(copy)
        tmp = tmp + AAs{i}(:,copy(j));
    end
    AAs{i}(:,orig) = AAs{i}(:,orig) + tmp;
end
