clear all

aX = 700;                                                                          
aH = 50;                                                                            
b = 0.000045:0.000017:0.0003;                                                                                                                            % Rate of Sup35 aggregation 
g = 1:0.5:13;                                                         

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
    dummy_aggrs = [1:98; repmat(20, [1,29]), repmat(40, [1,59]), repmat(50, [1,10])]';
    init_div_c = 1;

    while(init_div_c <= 4)

        [agg_r, ~] = size(dummy_aggrs);
        pp = rand([agg_r, 1]);
        init_matr = [];
        init_dummy = 1;

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

                init_matr = [init_matr; init_dummy, dummy_aggrs(xxx, 2)];
                init_dummy = init_dummy + 1;

            end

        end

        init_div_c = init_div_c + 1;
        dummy_aggrs = init_matr;

    end

    [agg_r, ~] = size(dummy_aggrs);
    pp = rand([agg_r, 1]);
    init_matr = [];
    init_dummy = 1;

    for(yyy = 1:agg_r)

        if(pp(yyy) >= 0.6)

            init_matr = [init_matr; init_dummy, dummy_aggrs(yyy, 2)];
            init_dummy = init_dummy + 1;

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

    for(j = 25:length(g))

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
            no = 0
            
            while(time_counter <= 1000)

                [agg_r, ~] = size(Y);

                if(agg_r == 0)

                    mycell(4) = mycell(4) + floor((mycell(2) - time_counter)*aX);
                    mycell(5) = mycell(5) + floor((mycell(2) - time_counter)*aH);
                    time_counter = mycell(2);
                    
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
                            
                            %HERE IS ONE BIG BIG
                            %PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                            if((cur_len - place) >= 6)

                                  Y(cur_agg, 2) = cur_len - place;

                            else

                                  mycell(4) = mycell(4) + cur_len - place;
                                  Y = [Y(1:

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

        fraction_plus = cur_psi_plus/k
        final_sup_ratio = mean(cur_ratio_sup)
        ratio_matrix(26-j, i) = fraction_plus;
        sup_ratio_matrix(26-j, i) = final_sup_ratio;

    end  

end

if(m == 1)

    our_case = 'no initial cell';

else

    our_case = 'with initial cell';

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
title(strcat('Prion stability, ', our_case))

figure
colormap('jet')
imagesc(sup_ratio_matrix)
caxis([0 1])
ax2 = gca;
ax2.XTick = [];
ax2.YTick = [];
xlabel('Conversion (b), ascending')
ylabel('Fragmentation(g), ascending')
title(strcat('Soluble Sup35 proportion, ', our_case))

