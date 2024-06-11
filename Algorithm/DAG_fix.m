function ind_vec = dag_fix(obj,ind_vec) 
            if any(ind_vec) == 0
                ind_vec = round(rand(1, size(ind_vec,2)));
            end
            interLength = max(ind_vec)-min(ind_vec); if interLength==0; interLength = 0.1;end
            ind_vec = round((ind_vec-min(ind_vec))/interLength);
            ind_adja = obj.vec_to_adja(ind_vec);
            ind_adja=obj.remove_cycles(ind_adja);
         
            ind_vec = obj.adja_to_vec(ind_adja);
            
end
        function DAG = remove_cycles(obj,adj_matrix)
            [rows, ~] = size(adj_matrix); 
            visited = zeros(rows, 1);
            on_stack = zeros(rows, 1); 
            cycle_detected = false;
            stack = [];
            for vertex = 1:rows
                if visited(vertex) == 0
                    visited(vertex) = 1;
                    on_stack(vertex) = 1;
                   stack = [stack vertex];
                    while ~isempty(stack)
                        current_vertex = stack(end);
                        adj_vertices = find(adj_matrix(current_vertex, :)); 
                        unvisited_adj_vertices = adj_vertices(visited(adj_vertices) == 0);
                        if ~isempty(unvisited_adj_vertices)
                            next_vertex = unvisited_adj_vertices(1);
                            visited(next_vertex) = 1; 
                            on_stack(next_vertex) = 1;
                            stack = [stack next_vertex]; 
                        else
                            stack = stack(1:end-1); 
                            on_stack(current_vertex) = 0; 
                            cycle_detected = any(on_stack(stack));
                            if cycle_detected
                                break;
                            end
                        end
                    end
                    if cycle_detected
                        break;
                    end
                end
            end

            if cycle_detected
                cycle = unique(stack); 
                for i = 1:numel(cycle)-1 
                    adj_matrix(cycle(i), cycle(i+1)) = 0; 
                end
                adj_matrix(cycle(end), cycle(1)) = 0; 
            end

            DAG = adj_matrix;
        end