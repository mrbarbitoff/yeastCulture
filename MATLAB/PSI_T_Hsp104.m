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

T = 20;											% Threshold level for aggregate inheritance

ratio_matrix = repmat(0, [length(g), length(b)]);
sup_ratio_matrix = repmat(0, [length(g), length(b)]);

for(j = 1:length(b))
    cur_b = b(j);
    for(k = 1: length(g))
        cur_g = g(k);
        
%         This piece of shitcode repeats five division cycles of the starting cell.
%         It repeats 4 cycles, evaluating the number of Sup35p copies (X), Hsp104 hex-s (H) and aggregate distribution (Y) (not generating 'em).
%         During this four we are considering only the start cell.
%         The fifth cycle generates the daughter cell, that is passed to 'cells' array with div time, X, H and Y{1}
        
        test_st = 1;

        while(test_st == 1)
            X = 30000;                                                                          
            H = 2000;
            Y = {};
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
            
            [test_r, test_c] = size(Y{1});
            
            if(test_r < 2)
                test_st = 1;
            else
                test_st = 0;
            end    
        end
        
        '[PSI+] cell with 2+ aggregates is successfully generated!'
        
%         May be useful to check for 2+ aggregates within a cell
        
%         HERE the initial cell for sim is finally generated and we are at last ready to begin it!
%         Each time step (i) the cell is checked for division (if it is dividing - it divides! XD);
%         the number of Sup35 and Hsp1-4 copies is computed, then if the cell is already [psi-] it is skipped,
%         if it is [PSI+] it undergoes several Gillespie cycles (until the time_counter reaches 1 minute).
        
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
                                if(Y{n}(o, 2) <= T)
                                    dummy_d = dummy_d + 1;
                                    matr_d = [matr_d; dummy_d, Y{n}(o, 2)];
                                else 
                                    dummy = dummy + 1;
                                    matr = [matr; dummy, Y{n}(o, 2)];
                                end    
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
                            prop_disagg = [];
                            for(p = 1:agg_r)
                                prop_agg = cells(n, 3)*cur_b;
                                prop_disagg = [prop_disagg; p, ((Y{n}(p, 2) - 1)*cur_g*cells(n, 4))/(cells(n, 4)/2 + (Y{n}(p, 2) - 1)), Y{n}(p, 2) - 1];
                            end    
                            total_prop = prop_agg + sum(prop_disagg(:, 2));
                            interval = exprnd(1/total_prop);
                            reaction = rand();
                            if(reaction >= prop_agg/total_prop)
                                ran_agg = unidrnd(sum(prop_disagg(:, 3)));
                                cur_agg = 0;
                                agg_def = 0;
                                while(agg_def < ran_agg)
                                    cur_agg = cur_agg + 1;
                                    agg_def = agg_def + prop_disagg(cur_agg, 3);
                                end
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
        
        
%         After 1000 minutes have passed, we calculate the proportion of [PSI+] cells and write it to a cell.
%         (For respective values of conversion and fragmentation rates).
        
        cond = cellfun(@isempty, Y);
        cond_ok = ~cond;
        filter = cond_ok == 1;
        norm = cond_ok(filter);
        no_plus = length(norm);
        [cell_r, cell_c] = size(cells);
        fr_plus = no_plus/cell_r
        ratio_matrix(26-k, j) = fr_plus;

	sup_test = unidrnd(cell_r, [100, 1]);
	sup_test_ratios = [];
	for(r = 1:length(sup_test)
		st_cell = sup_test(r);
		sup_test_cell = cells(st_cell, :);
		soluble_sup = sup_test_cell(3);
		sup_test_agg = Y{st_cell}(:, 2);
		insoluble_sup = sum(sup_test_agg);
		sup_test_ratios = [sup_test_ratios, soluble_sup/(soluble_sup + insoluble_sup)]
	end
	final_ratio = mean(sup_test_ratios)
	sup_ratio_matrix(26-k, j) = final_ratio;
    end
end   

%         Here we are finally ready to construct a heatmap using 'jet' colormap (like in the original paper).

colormap('jet')
imagesc(ratio_matrix)

imagesc(sup_ratio_matrix)