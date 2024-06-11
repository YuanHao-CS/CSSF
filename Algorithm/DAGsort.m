function Data = RandSort(Data)
Data = Data(randperm(size(Data, 2)));
end
function Population = PopSort(Population,optimization_direction)
PopSize = length(Population);
for i=1:PopSize-1
    for j=1:PopSize-i  
        if (optimization_direction==1)&&(Population(1,j).obj>Population(1,j+1).obj)
            tmp=Population(j);Population(j)=Population(j+1);Population(j+1)=tmp;
        end
        if (optimization_direction==-1)&&(Population(1,j).obj<Population(1,j+1).obj)
            tmp=Population(j);Population(j)=Population(j+1);Population(j+1)=tmp;
        end
    end
end
end