clear all

aX = 700;                                                                          
aH = 50;                                                                            
b = [0.000045:0.000017:0.0003];                                                                                                                            % Rate of Sup35 aggregation 
g = [0.0001:0.00005:0.0013];                                                         

X = 30000;
H = 2000;
aggrs = [1:98; repmat(20, [1,29]), repmat(40, [1,59]), repmat(50, [1,10])]';

rho = 31.03;                                                                            
lambda_d = 0.21*60;                                                                     
lambda = 1.16*60;                                                                       
shape_d = rho*lambda_d;                                                                 
shape = rho*lambda;                                                                      

ratio_matrix = repmat(0, [length(g), length(b)]);
sup_ratio_matrix = repmat(0, [length(g), length(b)]);

for(i = 1:length(b))
    
    cur_b = b(i);
    psi_plus = zeros(length(g), 1);
    ratio_sup = cell(length(g), 1)
    
    parfor(j = 1:length(g))
        
        cur_g = g(j);
        culture_Y = {};
        culture_Y{1} = aggrs;
        cells = [1, gamrnd(shape_d, 1/rho) + gamrnd(shape, 1/rho), 0, X, H];
        cur_psi_plus = psi_plus(j);
        cur_ratio_sup = ratio_sup{j};
                
        for(k = 1:2)
            
            mycell = cells(k, :);
            Y = culture_Y{k};
            time_counter = mycell(3);
            
            while(time_counter <= 1000)
                
                [agg_r, agg_c] = size(Y);
                
                if(agg_r == 0)
                    
                    time_counter = mycell(2);
                    
                end
                
                %Checking the division
                
                if(time_counter >= mycell(2))
                    
                    [cell_r, cell_c] = size(cells);
                    cells = [cells; cell_r + 1, mycell(2) + gamrnd(shape, 1/rho) + gamrnd(shape_d, 1/rho), mycell(2), 0.4*mycell(4), 0.4*mycell(5)];
                    mycell(2) = mycell(2) + gamrnd(shape, 1/rho);
                    mycell(4) = 0.6*mycell(4);
                    mycell(5) = 0.6*mycell(5);
                    time_counter = time_counter + 1;
                    
                    if(agg_r == 0)
                        
                        culture_Y{cell_r + 1} = [];
                        mycell(4) = (mycell(2) - time_counter)*aX;
                        mycell(5) = (mycell(2) - time_counter)*aH;
                        time_counter = mycell(2);
                        
                    else    
                        
                        pp = rand([agg_r, 1]);
                        matr = [];
                        matr_d = [];
                        dummy = 0;
                        dummy_d = 0;
                        
                        for(l = 1:agg_r)
                            
                            if(pp(l) >= 0.4)
                                
                                dummy = dummy + 1;
                                matr = [matr; dummy, Y(l, 2)];
                                
                            else
                                
                                dummy_d = dummy_d + 1;
                                matr_d = [matr_d; dummy_d, Y(l, 2)];
                                
                            end
                            
                        end
                        
                        Y = matr;
                        culture_Y{cell_r + 1} = matr_d;
                        
                    end
                    
                %If there is no division, behave as always    
                    
                else
                    
                    prop_agg = mycell(4)*agg_r*cur_b;
                    Z = sum(Y(:, 2)) - agg_r;
                    prop_disagg = (Z*cur_g*mycell(5))/((mycell(5)/2) + Z);
                    total_prop = prop_agg + prop_disagg;
                    interval = exprnd(1/total_prop);
                    mycell(4) = mycell(4) + aX*interval;
                    mycell(5) = mycell(5) + aH*interval;
                    reaction = rand();
                            
                    if(reaction >= prop_agg/total_prop)
                       
                        one_site_part = (prop_disagg/total_prop)/Z;
                        reaction_disagg = reaction - prop_agg/total_prop;
                        site = floor(reaction_disagg/one_site_part);
                        
                        if(site == 0)
                            
                            site = 1;
                            
                        end
                        
                        cur_agg = 0;
                        agg_def = 0;
                        
                        while(agg_def < site)
                            
                             cur_agg = cur_agg + 1;
                             agg_def = agg_def + Y(cur_agg, 2) - 1;
                             
                        end
                        
                        cur_len = Y(cur_agg, 2);
                        place = agg_def - site;
                            
                        if(place >= 6)
                            
                             Y(cur_agg, 2) = place;
                             
                             if((cur_len - place) >= 6)
                                 
                                  Y(agg_r + 1, :) = [agg_r + 1, cur_len - place];
                                  
                             else
                                 
                                  mycell(4) = mycell(4) + cur_len - place;
                                        
                             end
                             
                        else
                            
                            mycell(4) = mycell(4) + place;
                            
                            if((cur_len - place) >= 6)
                                
                                  Y(cur_agg, 2) = cur_len - place;
                                        
                            else
                                
                                mycell(4) = mycell(4) + cur_len - place;
                                        
                            end
                            
                        end     
                               
                    else
                        
                        cur_agg = unidrnd(agg_r);
                        Y(cur_agg, 2) = Y(cur_agg, 2) + 1;
                        mycell(4) = mycell(4) - 1;
                            
                    end
                    
                    time_counter = time_counter + interval;
                                  
                end
                               
            end
            
            cells(k, :) = mycell;
            culture_Y{k} = Y;
            psi_status = ~(isempty(Y));
            if(psi_status == 1)
                cur_psi_plus = cur_psi_plus + 1;
            end
            soluble_sup = mycell(4);
            insoluble_sup = sum(Y(:, 2));
            cur_ratio_sup = [cur_ratio_sup, soluble_sup/(soluble_sup + insoluble_sup)];
                      
        end
        
        fraction_plus = psi_plus/k
        final_sup_ratio = mean(ratio_sup{j})
        ratio_matrix(26-j, i) = fraction_plus;
        sup_ratio_matrix(26-j, i) = final_sup_ratio;
        
    end  
    
end    