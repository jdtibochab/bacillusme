%% Generate reaction ID dictionary
load('iYO844.mat')
rxnLink='http://bigg.ucsd.edu/api/v2/universal/reactions/';
filename = 'reaction_id_links.txt';
filename2 = 'reaction_id_notfound.txt';

%%
fileID = fopen(filename,'w');
fileID2 = fopen(filename2,'w');
model = iYO844;
notFound = {};
RXNS = model.rxns;

for i = 1:length(RXNS)
    rxnID = RXNS{i};
    rxnLink_rxn = [rxnLink,rxnID];
    try
        rxn(i)=webread(rxnLink_rxn);
        cycID = rxn(i).database_links.BioCyc.id;
        cycID = cycID(6:end);
        string = [rxnID, sprintf('\t'),cycID,sprintf('\n')];
        fprintf(fileID, string);
    catch
        cycID = 'NOT FOUND';
        string = [rxnID, sprintf('\t'),cycID,sprintf('\t'),model.rxnNames{i},sprintf('\t'),model.grRules{i},sprintf('\n')];
        fprintf(fileID2, string);
        notFound{end+1,1} = rxnID;
    end
end

%%
rxnIDs = findRxnIDs(model,dict);
grRules = model.grRules(rxnIDs);


