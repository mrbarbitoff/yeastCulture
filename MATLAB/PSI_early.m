%% Simulation with no Hsp104 and equal transmission of all aggregates

clear all

aX = 700;                                                                           % Rate of Sup35p biosynthesis mol/min
aH = 50;                                                                            % Rate of new Hsp104 hexamers formation mol/min
b = [0.000045:0.000017:0.0003];                                                     % Range of aggregation rates                                                                       % Rate of Sup35 aggregation 
g = [0.0001:0.00005:0.0013];                                                         % Rate of Hsp104-mediated disaggregation

% cur_b = 0.000045;
% cur_g = 0.0001;

rho = 31.03;                                                                            % Shape parameter of distribution, defines variance
lambda_d = 0.21*60;                                                                     % Mean time to become a "mother" cell
lambda = 1.16*60;                                                                       % Mean time to for a "mother" to produce a bud
shape_d = rho*lambda_d;                                                                 % Scale parameter of the "daughter-to-mother" distribution
shape = rho*lambda;                                                                     % Scale parameter of the 

ratio_matrix = repmat(0, [length(g), length(b)]);

for(j = 1:length(b))
    cur_b = b(j);
    for(k = 1: length(g))
        cur_g = g(k);
        X = 30000;                                                                          
        H = 2000;                                                                         
        aggrs = [1:98; repmat(20, [1,29]), repmat(40, [1,59]), repmat(50, [1,10])]';
        Y{1} = aggrs;
        a = 1;
        while(a <= 4)
            [agg_r, agg_c] = size(Y{1});
            pp = rand([agg_r, 1]);
            matr = [];
            dummy = 1;
            if (a == 1)
                tim = gamrnd(shape, 1/rho) + gamrnd(shape_d, 1/rho);
                X = 0.6*(X + tim*aX);
                H = 0.6*(H + tim*aH);
            else  
                tim = gamrnd(shape, 1/rho);
                X = 0.6*(X + tim*aX);
                H = 0.6*(H + tim*aH);
            end    
            for(l = 1:agg_r)
                if(pp(l) >= 0.4)
                    matr = [matr; dummy, Y{1}(l, 2)];
                    dummy = dummy + 1;
                end
            end
            a = a + 1;
            Y{1} = matr;
        end    
        [agg_r, agg_c] = size(Y{1});
        pp = rand([agg_r, 1]);
        matr = [];
        dummy = 1;
        for(m = 1:agg_r)
            if(pp(m) >= 0.6)
                matr = [matr; dummy, Y{1}(m, 2)];
                dummy = dummy + 1;
            end
        end
        Y{1} = matr;  
        tim = gamrnd(shape, 1/rho);
        X = 0.4*(X + tim*aX);
        H = 0.4*(H + tim*aH);
        cells = [1, gamrnd(shape, 1/rho) + gamrnd(shape_d, 1/rho), X, H];
        init = Y{1};
        
        % HERE the initial cell for sim is finally generated and we are finally ready to begin it!
        
        for(i = 1:1000)
            [cell_r, cell_c] = size(cells);
            ndiv = 0;
            for(n = 1:cell_r)
                time_counter = 0;
                mycell = cells(n, :);
                cells(n, 3) = cells(n, 3) + aX;
                cells(n, 4) = cells(n, 4) + aH;
                if(i >= mycell(2))
                    ndiv = ndiv + 1;
                    cells(n, 2) = mycell(2) + gamrnd(shape, 1/rho);
                    cells(n, 3) = 0.6*cells(n, 3);
                    cells(n, 4) = 0.6*cells(n, 4);
                    cells=[cells; cell_r + ndiv, mycell(2) + gamrnd(shape, 1/rho) + gamrnd(shape_d, 1/rho),0.4*cells(n, 3), 0.4*cells(n, 4)];
                    [agg_r, agg_c] = size(Y{n});
                    if(agg_r == 0)
                        Y{cell_r + ndiv} = [];
                    else    
                        pp = rand([agg_r, 1]);
                        matr = [];
                        matr_d = [];
                        dummy = 0;
                        dummy_d = 0;
                        for(o = 1:agg_r)
                            if(pp(o) >= 0.4)
                                dummy = dummy + 1;
                                matr = [matr; dummy, Y{n}(o, 2)];
                            else
                                dummy_d = dummy_d + 1;
                                matr_d = [matr_d; dummy_d, Y{n}(o, 2)];
                            end
                        end
                        Y{n} = matr;
                        Y{cell_r + ndiv} = matr_d;
                    end    
                else
                    while(time_counter <= 1)
                        [agg_r, agg_c] = size(Y{n});
                        if(agg_r == 0)
                            time_counter = 2;
                        else    
                            prop_agg = cells(n, 3)*agg_r*cur_b;
                            
%                             total_len = sum(Y{n}(:, 2));
%                             d_sites = total_len - agg_r;
%                             prop_disagg = d_sites*cur_g;
                            
                            prop_disagg = agg_r*cur_g;

                            total_prop = prop_agg + prop_disagg;
                            interval = exprnd(1/total_prop);
                            reaction = rand();
                            if(reaction >= prop_agg/total_prop)
                                cur_agg = unidrnd(agg_r);
                                cur_len = Y{n}(cur_agg, 2);
                                place = unidrnd(cur_len - 1);
                                if(place > 6)
                                    Y{n}(cur_agg, 2) = place;
                                    if((cur_len - place) >= 6)
                                        Y{n}(agg_r + 1, :) = [agg_r + 1, cur_len - place];
                                    else
                                        cells(n, 3) = cells(n, 3) + (cur_len - place);
                                    end
                                else
                                    cells(n, 3) = cells(n, 3) + place;
                                    if((cur_len - place) >= 6)
                                        Y{n}(cur_agg, 2) = cur_len - place;
                                    else
                                        cells(n, 3) = cells(n, 3) + (cur_len - place);
                                    end
                                end
                            else
                                cur_agg = unidrnd(agg_r);
                                Y{n}(cur_agg, 2) = Y{n}(cur_agg, 2) + 1;
                                cells(n, 3) = cells(n, 3) - 1;
                            end
                            time_counter = time_counter + interval;
                        end    
                    end    
                end  
            end   
        end
        cond = cellfun(@isempty, Y);
        cond_ok = ~cond;
        filter = cond_ok == 1;
        norm = cond_ok(filter);
        no_plus = length(norm);
        [cell_r, cell_c] = size(cells);
        fr_plus = no_plus/cell_r
        ratio_matrix(k, j) = fr_plus;
    end
end   

HeatMap(ratio_matrix, 'Colormap', 'redgreencmap')