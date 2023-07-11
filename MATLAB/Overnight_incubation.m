clear all

aX = 700;                                                                          
aH = 50;                                                                            
b = 0.0000045:0.0000197:0.0003;                                                                                                                            % Rate of Sup35 aggregation 
g = 0.001:0.0005:0.013;                                                         

rho = 31.03;                                                                            
lambda_d = 0.21*60;                                                                     
lambda = 1.16*60;                                                                       
shape_d = rho*lambda_d;                                                                 
shape = rho*lambda;                                                                                     

ratio_matrix = zeros([length(g), length(b)]);
sup_ratio_matrix = zeros([length(g), length(b)]);

test_st = 1;

while(test_st == 1)
    
    
    X = 30000;
    H = 2000;
    dummy_aggrs = [repmat(20, [1,29]), repmat(40, [1,59]), repmat(50, [1,10])];
    init_div_c = 1;

    while(init_div_c <= 4)

        agg_r = length(dummy_aggrs);
        pp = rand([agg_r, 1]);
        init_matr = [];

        if (init_div_c == 1)

            tim = gamrnd(shape, 1/rho) + gamrnd(shape_d, 1/rho);
            X = floor(0.6*(X + tim*aX));
            H = floor(0.6*(H + tim*aH));

        else  

            tim = gamrnd(shape, 1/rho);
            X = floor(0.6*(X + tim*aX));
            H = floor(0.6*(H + tim*aH));

        end  

        for(xxx = 1:agg_r)

            if(pp(xxx) >= 0.4)

                init_matr = [init_matr, dummy_aggrs(xxx)];

            end

        end

        init_div_c = init_div_c + 1;
        dummy_aggrs = init_matr;

    end

    agg_r = length(dummy_aggrs);
    pp = rand([agg_r, 1]);
    init_matr = [];

    for(yyy = 1:agg_r)

        if(pp(yyy) >= 0.6)

            init_matr = [init_matr, dummy_aggrs(yyy)];

        end

    end

    aggrs = init_matr;  
    tim = gamrnd(shape, 1/rho);
    X = floor(0.4*(X + tim*aX));
    H = floor(0.4*(H + tim*aH));

    test_st = isempty(aggrs);

end

for(i = 1:length(b))

    cur_b = b(i);
    psi_plus = zeros(length(g), 1);
    ratio_sup = cell(length(g), 1);

    parfor(j = 1:length(g))

        cur_g = g(j);
        culture_Y = {};
        culture_Y{1} = aggrs;
        cells = [1, gamrnd(shape_d, 1/rho) + gamrnd(shape, 1/rho), 0, X, H];
        cur_psi_plus = psi_plus(j);
        cur_ratio_sup = ratio_sup{j};

        for(k = 1:500)

            mycell = cells(k, :);
            Y = culture_Y{k};
            time_counter = mycell(3);
            
            while(time_counter <= 1000)

                agg_r = length(Y);

                if(agg_r == 0)

                    if(mycell(2) >= 1000)
                        
                        mycell(4) = mycell(4) + floor((1000 - time_counter)*aX);
                        mycell(5) = mycell(5) + floor((1000 - time_counter)*aH);
                        time_counter = 1000;
                        
                    else    
                        
                        mycell(4) = mycell(4) + floor((mycell(2) - time_counter)*aX);
                        mycell(5) = mycell(5) + floor((mycell(2) - time_counter)*aH);
                        time_counter = mycell(2);
                        
                    end
                    
                end

                %Checking the division

                if(time_counter >= mycell(2))

                    [cell_r, cell_c] = size(cells);
                    cells = [cells; cell_r + 1, mycell(2) + gamrnd(shape, 1/rho) + gamrnd(shape_d, 1/rho), mycell(2), floor(0.4*mycell(4)), floor(0.4*mycell(5))];
                    mycell(2) = mycell(2) + gamrnd(shape, 1/rho);
                    mycell(4) = floor(0.6*mycell(4));
                    mycell(5) = floor(0.6*mycell(5));
                    time_counter = time_counter + 1;

                    if(agg_r == 0)

                        culture_Y{cell_r + 1} = [];

                    else    

                        pp = rand([agg_r, 1]);
                        matr = [];
                        matr_d = [];

                        for(l = 1:agg_r)

                            if(pp(l) >= 0.4)

                                matr = [matr, Y(l)];
                                
                            else
                                
                                matr_d = [matr_d, Y(l)];

                            end

                        end

                        Y = matr;
                        culture_Y{cell_r + 1} = matr_d;

                    end

                %If there is no division, behave as always    

                elseif(time_counter == 1000)
                    
                    time_counter = 1001;
                                
                else

                    prop_agg = mycell(4)*agg_r*cur_b;
                    Z = sum(Y) - agg_r;
                    prop_disagg = (Z*cur_g*mycell(5))/((mycell(5)/2) + Z);
                    total_prop = prop_agg + prop_disagg + aX + aH;
                    interval = exprnd(1/total_prop);
                    reaction = rand();
                    
                    if(reaction <= aX/total_prop)
                        
                        mycell(4) = mycell(4) + 1;
                    
                    elseif(reaction <= (aX + aH)/total_prop)
                        
                        mycell(5) = mycell(5) + 1;
                      
                    elseif(reaction <= (aX + aH + prop_disagg)/total_prop)

                        one_site_part = (prop_disagg/total_prop)/Z;
                        reaction_disagg = reaction - (aX + aH)/total_prop;
                        site = floor(reaction_disagg/one_site_part);

                        if(site == 0)

                            site = 1;

                        end

                        cur_agg = 0;
                        agg_def = 0;

                        while(agg_def < site)

                             cur_agg = cur_agg + 1;
                             agg_def = agg_def + Y(cur_agg) - 1;

                        end

                        cur_len = Y(cur_agg);
                        place = agg_def - site;

                        if(place >= 6)

                             Y(cur_agg) = place;

                             if((cur_len - place) >= 6)

                                  Y(agg_r + 1) = cur_len - place;

                             else

                                  mycell(4) = mycell(4) + cur_len - place;

                             end

                        else
                            
                            mycell(4) = mycell(4) + place;

                            if((cur_len - place) >= 6)

                                  Y(cur_agg) = cur_len - place;

                            else

                                  mycell(4) = mycell(4) + cur_len - place;
                                  Y = [Y(1:(cur_agg - 1)), Y((cur_agg + 1):end)];

                            end

                        end     

                    else

                        cur_agg = unidrnd(agg_r);
                        Y(cur_agg) = Y(cur_agg) + 1;
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
            insoluble_sup = sum(Y);
            cur_ratio_sup = [cur_ratio_sup, soluble_sup/(soluble_sup + insoluble_sup)];

        end

        fraction_plus = cur_psi_plus/k
        final_sup_ratio = mean(cur_ratio_sup)
        ratio_matrix(26-j, i) = fraction_plus;
        sup_ratio_matrix(26-j, i) = final_sup_ratio;

    end  

end

figure
colormap('jet')
imagesc(ratio_matrix)
caxis([0 1])
ax1 = gca;
ax1.XTick = [];
ax1.YTick = [];
xlabel('Conversion (b), ascending')
ylabel('Fragmentation(g), ascending')
title('Prion stability, with initial cell')

figure
colormap('jet')
imagesc(sup_ratio_matrix)
caxis([0 1])
ax2 = gca;
ax2.XTick = [];
ax2.YTick = [];
xlabel('Conversion (b), ascending')
ylabel('Fragmentation(g), ascending')
title('Soluble Sup35 proportion, with initial cell')

