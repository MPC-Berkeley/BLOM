function AAs = MoveEqualVar(AAs,orig,copy)

for i=1:length(AAs)
    AAs{i}(:,orig) = AAs{i}(:,orig) + AAs{i}(:,copy);
end
