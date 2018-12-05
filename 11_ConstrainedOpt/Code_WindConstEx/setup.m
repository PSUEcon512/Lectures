% This is the main setup file
% Load data
function out = setup(varargin) 
%Input argument is a switch, when it is 0, we calculate road
% distances without allowing firms to use ferries, when it is 1 (the default), we
% calculate road distances using ferries (i.e., straight from google).
    if size(varargin)==0;
        use_ferries = 1;
    else
        use_ferries = varargin{1};
    end

    % old data using great cirlce distance
    % german_data = importdata('german_data_for_matlab.out',',',1);
    % danish_data = importdata('danish_data_for_matlab.out',',',1);

    % new data using road distances from man to producers and minimum great cirlce
    % distances to the 5 border crossing points
    german_data = importdata('german_data_for_matlab_new.out',',',1);
    if use_ferries
        danish_data = importdata('danish_data_for_matlab_new.out',',',1);
    else
        danish_data = importdata('danish_data_for_matlab_new_roadonly.out', ',', 1);
    end

    % There are 10 firms plus the fringe active in Germany, and 6 firms plus
    % the fringe in DK

    num_firms_ger = 10;
    num_firms_dnk = 6;

    %Number of foreign firms active in each country
    num_for_ger=5;
    num_for_dnk=1;

    num_pro_ger   = size(german_data.data(:,1),1);
    num_pro_dnk   = size(danish_data.data(:,1),1);

    m.y_matrix_ger  = german_data.data(:,4:14);
    % order of firms in GER:
    %'fringe','bonus','nordtank','micon','vestas','windworld','nordex','enercon
    %','fuhrlaender','suedwind','tacke' 

    m.y_matrix_dnk  = danish_data.data(:,4:10);
    % order of firms in DNK:
    % 'fringe','bonus','nordtank','micon','vestas','windworld','nordex'

    m.foreign_mat_ger = [ones(num_pro_ger,5),zeros(num_pro_ger,5)];
    m.foreign_mat_dnk = [zeros(num_pro_dnk,5),ones(num_pro_dnk,1)];

    %Scale distance to 100 kilometers...
    m.dist_mat_ger = log( (german_data.data(:,15:24))/100 );
    % order of firms
    %'distance_bonus','distance_nordtank','distance_micon','distance_vestas','distance_windworld',
    %'distance_nordex','distance_enercon','distance_fuhr','distance_suedwind','distance_tacke';

    %Scale distance to 100 kilometers...
    m.dist_mat_dnk = log( (danish_data.data(:,15:20))/100 );
    % order of firms
    % {'distance_bonus','distance_nordtank','distance_micon','distance_vestas',
    % 'distance_windworld','distance_nordex';}

    % Distance to border
    m.dist_bor_ger = (german_data.data(:,end))/100;
    m.dist_bor_dnk = (danish_data.data(:,end))/100;

    % Total MW of project

    m.ger_mw = german_data.data(:,2);    %num_pro_ger x 1
    m.dnk_mw = danish_data.data(:,2);    %num_pro_dnk x 1

    % Structure

    m.num_firms_ger = num_firms_ger;
    m.num_firms_dnk = num_firms_dnk;
    m.num_pro_ger   = num_pro_ger;
    m.num_pro_dnk   = num_pro_dnk;
    m.num_for_ger = num_for_ger;
    m.num_for_dnk = num_for_dnk;
    
    %% Add state-level information in Germany, this file comes from alex, old code simply doesn't use these variables.
   % german_firmstate_dummies = csvread('stateDumSorted.csv',1,1);
    %These are dummies if a product crosses state lines. They need to be 
    %ordered appropriately to fit with the columns of m. The first five
    %firms are danish, so all zeros, the remained are ordered
    %'distance_nordex','distance_enercon','distance_fuhr','distance_suedwind','distance_tacke';
    %which are columns 4 2 6 5 3 in the inported matrix. 
%     m.state_mat_ger= [zeros(m.num_pro_ger,5), ...
%                     german_firmstate_dummies(:,4), ...
%                     german_firmstate_dummies(:,2), ...
%                     german_firmstate_dummies(:,6), ...
%                     german_firmstate_dummies(:,5), ...
%                     german_firmstate_dummies(:,3) ];   
    
    m.state_mat_ger= [zeros(m.num_pro_ger,5), german_data.data(:,end-5:end-1)];
    %This is not used, since there are no internal Danish borders in the
    %models
    m.state_mat_dnk = zeros(m.num_pro_dnk, m.num_firms_dnk);


    %% Create structure for no fixed cost counterfactuals...similar to the m structure but with every firm competing 
    %  in both markets..
    %Add the german domestic-only firms to the danish market
    m_nFC = m;
    m_nFC.num_firms_dnk = m.num_firms_ger;
    m_nFC.dist_mat_dnk = log((danish_data.data(:,15:24))/100);
    m_nFC.num_for_dnk = 5;
    m_nFC.foreign_mat_dnk = [zeros(num_pro_dnk,5),ones(num_pro_dnk,5)];
    m_nFC.state_mat_dnk = zeros(m_nFC.num_pro_dnk, m_nFC.num_firms_dnk); %There are no states in Denmark.
    
    %This is the vector of ouput variables e_% are the estimated values %_mat
    %are the bootstrap runs for each one...
    name_vec = {'rhohat' , 'rho_nF' ,  'rho_nB' , 'natShare', 'natShareNF' , 'natShareNB' , 'markups' ,  'markupsNF' , ...
        'markupsNB' , 'profits' ,  'profitsNF' , 'profitsNB' ,  'conSurp' , 'conSurpNF' , 'conSurpNB' , ...
        'delasticity_ger' , 'delasticity_dnk' , 'fc_bounds' };

    save('data_structures', 'm', 'm_nFC', 'name_vec', 'use_ferries');

    out = use_ferries;

end