clear all
close all
clc

A = importdata('NC_045512_site_database_altA.tsv'); %col 18:end are intrahost samples
A = cellfun(@(s) strsplit(s), A, 'UniformOutput', false);

T = importdata('NC_045512_site_database_altT.tsv');
T = cellfun(@(s) strsplit(s), T, 'UniformOutput', false);

G = importdata('NC_045512_site_database_altG.tsv');
G = cellfun(@(s) strsplit(s), G, 'UniformOutput', false);

C = importdata('NC_045512_site_database_altC.tsv');
C = cellfun(@(s) strsplit(s), C, 'UniformOutput', false);


row = size(A,1);
Result = zeros(row-1,8); %in the order of ATGC.
All = zeros(row-1,2);
ListNuc = [{'A'},{'T'},{'G'},{'C'}];
Gene1 = cellfun(@(s) s{3}, A(2:end), 'UniformOutput',false);
Gene2 = cellfun(@(s) s{4}, A(2:end), 'UniformOutput',false);
uniqueGene1 = setdiff(unique(Gene1),'.');
uniqueGene2 = setdiff(unique(Gene2),'.');

Coor1 = getGeneCoor(Gene1,uniqueGene1);
Coor2 = getGeneCoor(Gene2,uniqueGene2);


for i = 2:row
    recA= returnRecurrent(A{i}(18:end));
    recT = returnRecurrent(T{i}(18:end));
    recG = returnRecurrent(G{i}(18:end));
    recC = returnRecurrent(C{i}(18:end));
    Result(i-1,:) = [recA,recT,recG,recC];
    Data = [A{i}(18:end);T{i}(18:end);G{i}(18:end);C{i}(18:end)];
    indNuc = find(strcmp(ListNuc,A{i}(2)) == 1);
    Data(indNuc,:) = [];
    All(i-1,:) = returnRecurrent(reshape(Data,1,401*3));
end
indA = find(Result(:,1) >= 10);
indU = find(Result(:,3) >= 10);
nonsynAll = [{getNonSyn(A)},{getNonSyn(T)},{getNonSyn(G)},{getNonSyn(C)}];

Colors = [    0.5 0.5 0;
    0 0.5 0.5;
    0.5 0 0.5;
    0.1 0.1 0.1;
    0.1 0.5 0.7;
    0.3 0.7 0.5;
    1 0 0;
    0 0 1;
    0 1 0;

    1 1 0;
    0 1 1;
    1 0 1];

%-----plot the distribution----
subplot('Position',[0.08 0.2 0.86 0.75])
hold on

X = All([1:row-1]',1)./All([1:row-1]',2);
indX = find(X>0.025);
scatter(indX,X(indX),130,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','k');


freqExceed = [];
indExceed = [];
for i = 1:4;
    X = Result([1:row-1]',i*2-1)./All([1:row-1]',2);
    indX = find(X>0.025);
    scatter(indX,X(indX),130,'filled' ,'MarkerFaceAlpha',0.5); 
    indNS = nonsynAll{i};
    indExceed = [indExceed;intersect(indX, indNS)];
    freqExceed = [freqExceed;X(intersect(indX, indNS))];
end
scatter(indExceed, freqExceed,20,'filled','MarkerFaceColor','r');
ylim([0 .9])


indAll = find(All(:,1)./All(:,2) >= 0.025);
length(indAll);

[freqSort indSort] = sort(freqExceed);
indSortTop = indExceed(indSort(end-4:end))
freqSortTop = freqSort(end-4:end);
coorChange = [{'C10083U'},{'C5768A'},{'C8655U'},{'A14553C'},{'A404U'}]
aaChange = [{'S>F'},{'H>N'},{'S>F'},{'L>F'},{'K>Stop'}]
for i = 2:length(indSortTop)
    text(indSortTop(i)-200, freqSortTop(i) + 0.05,strcat(coorChange{i},', ',aaChange{i}),'FontSize',18);
end
text(10087-200,0.192,'T10087G','FontSize',18);
for i = 1:length(uniqueGene1)
    for j = 0:10:(Coor1(i,2) - Coor1(i,1))
    plot([Coor1(i,1)+j,Coor1(i,1)+j],[0.0,0.02],'Color',Colors(i,:),'HandleVisibility','off')
    end
end

for i = 1:length(uniqueGene2)
    for j = 0:5:(Coor2(i,2) - Coor2(i,1))
    plot([Coor2(i,1)+j,Coor2(i,1)+j],[0.005,0.02],'k','HandleVisibility','off')
    end
end

legend('All','A','U','G','C','Nonsyn');
xlabel('Genome position')
ylabel('Proportion of samples');
set(gca,'FontSize',22);
box on
set(gcf,'PaperPosition',[0 0 16 6]);
saveas(1,'Fig8.jpg');



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
