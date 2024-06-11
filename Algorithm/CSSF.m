classdef CSSF< ALGORITHM
% <single> <real> <large/none> <expensive/none>
% CEO: Co-Evolutionary Optimization
% ambient_pressure   --- 0.4 --- Ambient Pressure between [0,1]
% optimization_direction   --- -1 --- Use -1 and 1 respectively to represent the optimization direction of minimization or maximization

%------------------------------- Reference --------------------------------
% Tian-Yi Liu, Yuan-Hao Jiang, Yuang Wei, et al., Recognition Real-World Educational Learning Path with Collaborative Structural Search Framework, 2024
%------------------------------- Copyright --------------------------------
% The proposed algorithm CSSF is an add-on in the PlatEMO platform. Therefore, users are requested to follow the requirements of the PlatEMO platform. More information about the PlatEMO platform can be found here:

%Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform for evolutionary multi-objective optimization [educational forum], IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87
%----------------------------------------------------------------------   ----

% This function is written by Tianyi Liu,
% E-mail: 231210704104@stu.just.edu.cn

     methods
        function main(Algorithm,Problem)
            [ambient_pressure,optimization_direction] = Algorithm.ParameterSet(0.4,-1);
            
            Population = Problem.Initialization();
            Archive    = [];
            count_ones = zeros(1, length(Population(1,1).dec));
            MCR = zeros(Problem.N,1) + 0.5;
            MF  = zeros(Problem.N,1) + 0.5;
            k   = 1;
            PF=0.6;
            NF=0.4;
            N_min=5;
            N_init=Problem.N;
    
            terminal_value=100000;
            while Algorithm.NotTerminated(Population)
    
                Population = PopSort(Population,optimization_direction);
                PopSize=size(Population,2);
                [~,rank] = sort(FitnessSingle(Population));
                Xpb = Population(rank(ceil(rand(1,PopSize).*max(2,rand(1,PopSize)*0.2*PopSize)))).decs;
                Xr1 = Population(randi(end,1,PopSize)).decs;
                P   = [Population,Archive];
                Xr2 = P(randi(end,1,PopSize)).decs;
                CR  = randn(PopSize,1).*sqrt(0.1) + MCR(randi(end,PopSize,1));
                CR  = repmat(max(0,min(1,CR)),1,Problem.D);
                F   = min(1,trnd(1,PopSize,1).*sqrt(0.1) + MF(randi(end,PopSize,1)));
                while any(F<=0)
                    F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
                end
                F = repmat(F,1,Problem.D);
                Site         = rand(size(CR)) < CR;
                OffDec       = Population.decs;
                OffDec(Site) = OffDec(Site) + F(Site).*(Xpb(Site)-OffDec(Site)+Xr1(Site)-Xr2(Site));
                Offspring    = Problem.Evaluation(OffDec);
                Offspring=PopSort(Offspring,optimization_direction);


               for i=ambient_pressure*PopSize:PopSize 
                   NewEdge=Offspring(1,i).dec-Population(1,i).dec;
                   NewEdge=NewEdge>0;
                   count_ones=NewEdge+count_ones;
               end
               nonzero_idx = count_ones ~= 0;
               count_ones(nonzero_idx) = (count_ones(nonzero_idx) - min(count_ones(nonzero_idx))) / (max(count_ones(nonzero_idx)) - min(count_ones(nonzero_idx)));
               count_ones=count_ones*PF;
               for i=ambient_pressure*PopSize:PopSize
                   randnum=rand(1,length(Offspring(1,i).dec));
                   replace=count_ones-randnum>0;
                   tmp=Offspring(1,i).dec;
                   tmp(replace)=1;
                   Offspring(1,i)=Problem.Evaluation(tmp);
               end
            
               for i=1:ambient_pressure*PopSize-1
                   NewEdge=Offspring(1,i).dec-Population(1,i).dec;
                   NewEdge=NewEdge>0;
                   count_ones=NewEdge+count_ones;
               end
               nonzero_idx = count_ones ~= 0;
    
               count_ones(nonzero_idx) = (count_ones(nonzero_idx) - min(count_ones(nonzero_idx))) / (max(count_ones(nonzero_idx)) - min(count_ones(nonzero_idx)));
               count_ones=count_ones*NF;
               for i=1:ambient_pressure*PopSize-1
                   randnum=rand(1,length(Offspring(1,i).dec));
                   replace=count_ones-randnum>0;
                   tmp=Offspring(1,i).dec;
                   tmp(replace)=0;
                   Offspring(1,i)=Problem.Evaluation(tmp);
               end
      
                Offspring=PopSort(Offspring,-optimization_direction);
                delta   = FitnessSingle(Population) - FitnessSingle(Offspring);
                replace = delta > 0;
                Archive = [Archive,Population(replace)];
                Archive = Archive(randperm(end,min(end,PopSize)));
                Population(replace) = Offspring(replace);  
                if any(replace)
                    w      = delta(replace)./sum(delta(replace));
                    if MCR(k)==terminal_value||isempty(CR(replace))
                        MCR(k)=terminal_value;
                    else
                        MCR(k)=(w'*CR(replace,1).^2)./(w'*CR(replace,1));
                    end
                    MF(k)=(w'*F(replace,1).^2)./(w'*F(replace,1));
                    k = mod(k,length(MCR))+1;
    
                end
                Population = PopSort(Population,optimization_direction);
                Population(1:ambient_pressure*PopSize) = Population((PopSize-ambient_pressure*PopSize+1):PopSize);
    
                N_G=max(0,round((N_min-N_init)/(Problem.maxFE)*Problem.FE)+N_init);
                if N_G<size(Population,2)
                    Archive=Archive(randperm(end,min(end,N_G)));
                end
    
            end
        end
    end
end

