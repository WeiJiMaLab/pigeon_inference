subs = {'1_MG' '2_YZ' '3_CH' '4_NJ' '5_IK' '6_OM' '7_JL' '8_CM' '9_BB' '10_JM'};

for subjidx = 1:length(subs)
    % Compiles blocks for an individual subject and saves in folder
    % z_concat_data_expt1
    sub_file_ID = subs{subjidx};

    %% compile blocks
    addpath('/Users/jennlauralee/Google Drive/WedMock/Causal Inference/Pigeon_expt1/Functions')
    path_in = ['/Users/jennlauralee/Google Drive/WedMock/Causal Inference/Analysis/z_raw_data_expt1/' sub_file_ID];
    path_out = ['/Users/jennlauralee/Google Drive/WedMock/Causal Inference/Analysis/z_concat_data_expt1/'];

    cd(path_in)

    scrape  = dir('*.mat');
    nBlocks = size(scrape);

    Data_ = [];
    Stim_ = [];
    load(scrape(1).name);

    datafields = fieldnames(data);
    stimfields = fieldnames(stim);

    for i_block = 1:nBlocks
        load(scrape(i_block).name)
        Data_ = MergeStructs(Data_, data);
        Stim_ = MergeStructs(Stim_, stim);
    end

    Stim.Feeder =[];
    Stim.X = {};
    Stim.Y = {};
    Stim.N = [];
    Stim.N1 = [];
    Stim.N0 = [];
    Stim.SX = [];
    Stim.SY = [];
    Stim.X1 = {};
    Stim.Y1 = {};

    Data.Resp = [];
    Data.Resp_feeder = [];
    Data.Resp_conf = [];
    Data.Correct_feeder = [];
    Data.RT = [];

    for i_block = 1:nBlocks
        Stim.Feeder = [Stim.Feeder; Stim_(i_block).feeder];
        Stim.X = [Stim.X Stim_(i_block).x];
        Stim.Y = [Stim.Y Stim_(i_block).y];
        Stim.N = [Stim.N; Stim_(i_block).n];
        Stim.N1 = [Stim.N1; Stim_(i_block).n1];
        Stim.N0 = [Stim.N0; Stim_(i_block).n0];
        Stim.SX = [Stim.SX Stim_(i_block).sx];
        Stim.SY = [Stim.SY Stim_(i_block).sy];
        Stim.X1 = [Stim.X1 Stim_(i_block).x1];
        Stim.Y1 = [Stim.Y1 Stim_(i_block).y1];

        Data.Resp = [Data.Resp; Data_(i_block).resp'];
        Data.Resp_feeder = [Data.Resp_feeder; Data_(i_block).resp_feeder'];
        Data.Resp_conf = [Data.Resp_conf; Data_(i_block).resp_conf'];
        Data.Correct_feeder = [Data.Correct_feeder; Data_(i_block).correct_feeder'];
        Data.RT = [Data.RT; Data_(i_block).r1_RT'];
    end

    % Stimulus post-processing
    Stim.MeanDist = getMeanDist(Stim);

    save([path_out sub_file_ID '.mat'], 'Data', 'Stim');

    % Superstructure alldata
    DATA{subjidx} = Data;
    STIM{subjidx} = Stim;
end

save('/Users/jennlauralee/Google Drive/WedMock/Causal Inference/Analysis/alldata.mat', 'DATA', 'STIM')