function [correl_mat] = correl_calc(H,HPSflag)
%This section calculates the correlation coefficients (R^2) of the Hamiltonians and the functional data.  It also calculates the percentile of the R^2.

load('~/Documents/Postdoc/Faruck/Hana/get_alignments/results_tr_sp_TEAS_HPS_based_seed_final/file5_HPSData.mat');
load('~/Documents/Postdoc/Faruck/Hana/get_alignments/results_tr_sp_TEAS_HPS_based_seed_final/file6_TEASData.mat');

if(HPSflag == 1)
    [R,P] = r_calc(H,HPS_Prem);
    correl_mat(1,:) = {'Premnaspirodiene', round(R,3), round(P,2)};
    [R,P] = r_calc(H,HPS_5EA);
    correl_mat(2,:) = {'5-epi-aristolochene', round(R,3), round(P,2)};
    [R,P] = r_calc(H,HPS_4EE);
    correl_mat(3,:) = {'4-epi-eremophiene', round(R,3), round(P,2)};
    [R,P] = r_calc(H,HPS_GermA);
    correl_mat(4,:) = {'Germacrene A', round(R,3), round(P,2)};
    [R,P] = r_calc(H,HPS_Sol);
    correl_mat(5,:) = {'Solubility', round(R,3), round(P,2)};
    [R,P] = r_calc(H,HPS_Tm);
    correl_mat(6,:) = {'Tm', round(R,3), round(P,2)};
    [R,P] = r_calc(H,HPS_entropy);
    correl_mat(7,:) = {'Entropy', round(R,3), round(P,2)};
else
    [R,P] = r_calc(H,TEAS_Prem);
    correl_mat(1,:) = {'Premnaspirodiene', round(R,3), round(P,2)};
    [R,P] = r_calc(H,TEAS_5EA);
    correl_mat(2,:) = {'5-epi-aristolochene', round(R,3), round(P,2)};
    [R,P] = r_calc(H,TEAS_4EE);
    correl_mat(3,:) = {'4-epi-eremophiene', round(R,3), round(P,2)};
    [R,P] = r_calc(H,TEAS_GermA);
    correl_mat(4,:) = {'Germacrene A', round(R,3), round(P,2)};
    [R,P] = r_calc(H,TEAS_Sol);
    correl_mat(5,:) = {'Solubility', round(R,3), round(P,2)};
    [R,P] = r_calc(H,TEAS_Tm);
    correl_mat(6,:) = {'Tm', round(R,3), round(P,2)};
    [R,P] = r_calc(H,TEAS_entropy);
    correl_mat(7,:) = {'Entropy', round(R,3), round(P,2)};
end
end

function [R,percentile] = r_calc(H_mutations,func_data)
H_mutations = H_mutations(~isnan(func_data));
func_data = func_data(~isnan(func_data));
func_data = func_data(~isnan(H_mutations));
H_mutations = H_mutations(~isnan(H_mutations));
H_mutations = H_mutations(~isinf(func_data));
func_data = func_data(~isinf(func_data));

R = corrcoef(H_mutations,func_data);
R = R(2);
R2 = R.^2;
[rand_R] = random_seq_correlation(H_mutations,func_data);
sorted_rand_R = sort(rand_R);
sorted_rand_R2 = sort(rand_R.^2);

num_below = 0;

for i = 1:length(sorted_rand_R2)
    if(sorted_rand_R2(i) < R2)
        num_below = num_below+1;
%     else
%         disp(sorted_rand_R2(i));
     end
end

percentile = (1-(length(sorted_rand_R2) - num_below)/length(sorted_rand_R2))*100;
end
