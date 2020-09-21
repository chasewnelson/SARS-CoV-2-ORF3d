
clear all
close all
clc

A = importdata('NC_045512_site_database_altA.tsv'); % 18:end are intrahost samples
A = cellfun(@(s) strsplit(s), A, 'UniformOutput', false);

T = importdata('NC_045512_site_database_altT.tsv');
T = cellfun(@(s) strsplit(s), T, 'UniformOutput', false);

G = importdata('NC_045512_site_database_altG.tsv');
G = cellfun(@(s) strsplit(s), G, 'UniformOutput', false);

C = importdata('NC_045512_site_database_altC.tsv');
C = cellfun(@(s) strsplit(s), C, 'UniformOutput', false);

orf3clof = 14408;
recurrFreqLOF = returnFreq(T{orf3clof+1}(18:end));
%
coorChange = [{'A404U'},{'A14553C'},{'C8655U'}, {'C5768A'}];
xbins = [0:0.01:0.5];

recurrFreq404 = returnFreq(T{404+1}(18:end)); %add 1 for header
reccurFreq14553 = returnFreq(C{14553+1}(18:end));
reccurFreq8655 = returnFreq(T{8655+1}(18:end));
reccurFreq5768 = returnFreq(A{5768+1}(18:end));

Result = [];
subplot(1,4,1)
h1 = histogram(recurrFreq404,xbins);
Result = [Result,h1.Values'];
ylabel('Occurrence');
set(gca,'FontSize',22);
title(coorChange{1})

xlabel('DAF within-host')
box on
subplot(1,4,2)
h2 = histogram(reccurFreq14553,xbins);
Result = [Result,h2.Values'];
set(gca,'FontSize',22);
xlabel('DAF within-host')
box on
title(coorChange{2})
subplot(1,4,3)
h3 = histogram(reccurFreq8655,xbins);
Result = [Result,h3.Values'];
set(gca,'FontSize',22);
box on
title(coorChange{3})
xlabel('DAF within-host')
subplot(1,4,4)
h4 = histogram(reccurFreq5768,xbins);
Result = [Result,h4.Values'];
xlabel('DAF within-host')
set(gca,'FontSize',22);
box on
title(coorChange{4})

set(gcf,'PaperPosition',[0 0 18 5]);
saveas(1,'Fig8_supp1.jpg');


function nonsynIndex = getNonSyn(List)    
    aa1 = cellfun(@(s) s{11}, List(2:end), 'UniformOutput',false);
    aa2 = cellfun(@(s) s{16}, List(2:end), 'UniformOutput',false);
    aa3 = cellfun(@(s) s{12}, List(2:end), 'UniformOutput',false);
    aa4 = cellfun(@(s) s{17}, List(2:end), 'UniformOutput',false);
    diffAA1 = find(strcmp(aa1,aa2)~=1);
    diffAA2 = find(strcmp(aa3,aa4)~=1);
    nonsynIndex = union(diffAA1,diffAA2);
end

function Coor = getGeneCoor(GeneList,GeneName);
    Coor = zeros(length(GeneName),2);
    for i = 1:length(GeneName);
        ind = find(strcmp(GeneList,GeneName{i})==1);
        Coor(i,:) = [min(ind),max(ind)];
    end
end

function recurrFreq = returnFreq(A)
    recurrFreq = zeros(length(A),1);
    indKeep = find(contains(A,',') == 1);
    if length(indKeep) > 0
        data = cellfun(@(s) strsplit(s, ','), A(indKeep),'UniformOutput',false);
        refAllele = cell2mat(cellfun(@(s) str2num(s{1}), data,'UniformOutput',false));
        altAllele = cell2mat(cellfun(@(s) str2num(s{2}), data,'UniformOutput',false));
        recurrFreq(indKeep) = altAllele./(refAllele + altAllele);
    end
end

function recurrVec = returnRecurrent(A)
    indKeep = find(contains(A,',') == 1);
    if length(indKeep) > 0
        data = cellfun(@(s) strsplit(s, ','), A(indKeep),'UniformOutput',false);
        refAllele = cell2mat(cellfun(@(s) str2num(s{1}), data,'UniformOutput',false));
        altAllele = cell2mat(cellfun(@(s) str2num(s{2}), data,'UniformOutput',false));
        alleleFreq = altAllele./(refAllele + altAllele);
        ind1 = intersect(find(alleleFreq < 0.5),find(alleleFreq > 0 ));
        ind2 = intersect(find(alleleFreq > 0.5),find(alleleFreq < 1 ));
        recurrVec = [length(ind1), 401 - length(indKeep)]; 
    else
        recurrVec = [0, 0];
    end
end
