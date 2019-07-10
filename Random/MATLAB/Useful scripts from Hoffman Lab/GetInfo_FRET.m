function SaveParams = GetInfo_FRET(folder)

% A program allowing for the manual input of all parameters
% for a FRET experiment.

if not(exist(fullfile(folder,['SaveParams_' folder '.mat']),'file')) % Manually input and save parameters used in the analysis
    
    SaveParams.folder = folder;
    SaveParams.num_exp = input('How many experimental groups do you have? ');
    SaveParams.exp_cell = cell(1,SaveParams.num_exp);
    for i = 1:SaveParams.num_exp
        SaveParams.exp_cell{i} = input('Enter an experimental group name \n(Ex. VinTS): ','s');
    end
    SaveParams.num_channel = 3;
    SaveParams.mag = input('What magnification were your images taken at (40x or 60x)? ','s');
    SaveParams.temperature = input('What temperature were your images taken at ("23C" or "37C")? ','s');
    SaveParams.ND = input('What ND filter was in when you took images ("ND100" or "ND50" etc...)? ','s');
    SaveParams.Achannel = input('What is your acceptor channel? ','s');
    SaveParams.FRETchannel = input('What is your FRET channel? ','s');
    SaveParams.Dchannel = input('What is your donor channel? ', 's');
    
    % Only for FRET
    SaveParams.bt = input('Calculate bleedthroughs (y or n)? ','s');
    if strcmpi(SaveParams.temperature,'23C')
        if strcmpi(SaveParams.mag,'60x')
            SaveParams.G = 2.65;
            SaveParams.k = 0.68;
        elseif strcmpi(SaveParams.mag,'40x')
            SaveParams.G = 3.08;
            SaveParams.k = 0.855;
        end
    elseif strcmpi(SaveParams.temperature,'37C')
        SaveParams.G = 2.95;
        SaveParams.k = 0.66;
    end
    if strcmpi(SaveParams.bt,'y');
        SaveParams.donor_pre = input('Enter donor image names (Ex. Teal): ','s');
        SaveParams.acceptor_pre = input('Enter donor image names (Ex. Venus): ','s');
        SaveParams.dthres = input('Only calculate BTs above donor intensity (~500): ');
        SaveParams.athres = input('Only calculate BTs above acceptor intensity (~800): ');
    elseif strcmpi(SaveParams.bt,'n');
        SaveParams.dthres = 500;
        SaveParams.athres = 800;
        if strcmpi(SaveParams.mag,'40x') && strcmpi(SaveParams.temperature,'23C') && strcmpi(SaveParams.ND,'ND50')
            SaveParams.abt = 0.26;
            SaveParams.dbt = 0.95;
        elseif strcmpi(SaveParams.mag,'40x') && strcmpi(SaveParams.temperature,'23C') && strcmpi(SaveParams.ND,'ND100')
            SaveParams.abt = 0.27;
            SaveParams.dbt = 0.95;
        elseif strcmpi(SaveParams.mag,'60x') && strcmpi(SaveParams.temperature,'23C') && strcmpi(SaveParams.ND,'ND50')
            SaveParams.abt = 0.24;
            SaveParams.dbt = 0.94;
        elseif strcmpi(SaveParams.mag,'60x') && strcmpi(SaveParams.temperature,'23C') && strcmpi(SaveParams.ND,'ND100')
            SaveParams.abt = 0.27;
            SaveParams.dbt = 0.96;
        elseif strcmpi(SaveParams.mag,'60x') && strcmpi(SaveParams.temperature,'37C') && strcmpi(SaveParams.ND,'ND50')
            SaveParams.abt = 0.29;
            SaveParams.dbt = 1.06;
        elseif strcmpi(SaveParams.mag,'60x') && strcmpi(SaveParams.temperature,'37C') && strcmpi(SaveParams.ND,'ND100')
            SaveParams.abt = 0.31;
            SaveParams.dbt = 1.10;
        end
    end
    SaveParams.correct = input('FRET correct (y or n)? ','s');
    if strcmpi(SaveParams.correct,'y')
        SaveParams.venus_thres = input('Set all pixels to zero below venus threshold (~100): ');
    end
    
    save(fullfile(folder,['SaveParams_' folder '.mat']),'-struct','SaveParams');
    
else % Load the parameter file and save variables as the different parts of it
    SaveParams = load(fullfile(folder,['SaveParams_' folder '.mat']));
end
