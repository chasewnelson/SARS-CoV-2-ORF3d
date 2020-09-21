clear all
close all
clc

[Num Text] = xlsread('SARSCOV2_MAFFT_processed_variants_MAF0.02_ManualCheck_Fig7.xlsx');
Time = Num(2:end,1:3);
Coor = Num(1,5:end)+4;
Location = Text(2:end,3);
ID = Text(2:end,1);
Genotype = Text(2:end,8:end);

indTimePass = intersect(intersect(find(Time(:,2)>0),find(Time(:,1)>0)),find(Time(:,3)>=2019)); 
indTimeExclude = setdiff([1:length(ID)],indTimePass);
GenotypeExcluded = Genotype(indTimeExclude,:);
TimeExcluded = Time(indTimeExclude,:);
LocationExcluded = Location(indTimeExclude,:);
Time = Time(indTimePass,:);

Location = Location(indTimePass);
ID = ID(indTimePass);
Genotype = Genotype(indTimePass,:);
DateNumber = datenum(Time(:,3),Time(:,1),Time(:,2));
DateNumber = DateNumber - min(DateNumber);
[DateNumber indSort] = sort(DateNumber);
Time = Time(indSort,:);
Location = Location(indSort);
ID = ID(indSort);
Genotype = Genotype(indSort,:);
lettersToExclude = [{'-'},{'K'},{'N'},{'R'},{'S'},{'Y'}];
lettersAllowed = [{'A'},{'T'},{'G'},{'C'}]

TimeString = strcat(num2str(Time(:,1)),' - ',num2str(Time(:,2)));
MonthDay = [];
for i = 1:length(TimeString)
    MonthDay = [MonthDay;{TimeString(i,:)}];
end

Ancestral = Genotype(1,:);
Derived = [];
indexDerived = cell(length(Coor),1);
combined = [];
bias = [];
DAF = [];
DateDAFirstOccurance = [];
TimeDAF = [];
LocationDAFirstOccurance = [];
DAFIndices = [];
indDAFAll = [];
DaysDAF = [];
for i = 1:length(Coor)
    letters = intersect(unique(Genotype(:,i)),lettersAllowed);
    indDerived = find(contains(letters, Ancestral{i}) ~= 1);
    Derived = [Derived, letters(indDerived)];
    indexDerived{i} = find(strcmp(Genotype(:,i),letters(indDerived)) == 1);
    combined = [combined, {strcat(Ancestral{i},num2str(Coor(i)), Derived{i})}];
    indDA = find(strcmp(Genotype(:,i),letters(indDerived))==1);
    indAA = find(strcmp(Genotype(:,i),Ancestral{i})==1);
    DAF = [DAF, length(indDA)/(length(indDA)+length(indAA))];
    bias = [bias,{strcat(Ancestral{i},Derived{i})}];
    DateDAFirstOccurance = [DateDAFirstOccurance; DateNumber(indDA(1:5))'];
    TimeDAF = [TimeDAF;MonthDay(indDA(1:5),:)'];
    LocationDAFirstOccurance = [LocationDAFirstOccurance; Location(indDA(1:5))'];
    DAFIndices = [DAFIndices;indDA(1:5)'];
    indDAFAll = [indDAFAll;{indDA}];
    DaysDAF = [DaysDAF;DateNumber(indDA(1:5))'];
end

indGerman = [1,8,11,19,20];
indExcludeLetter = [];
GenotypeGerman = Genotype(:,indGerman);
for i = 1:length(indGerman);
    for j = 1:length(lettersToExclude)
        indExcludeLetter = [indExcludeLetter; find(strcmp(GenotypeGerman(:,i),lettersToExclude{j})==1)];
    end
end

GenotypeGerman(indExcludeLetter,:) = [];
MonthDay(indExcludeLetter) = [];
DateNumber(indExcludeLetter) = [];
Time(indExcludeLetter,:) = [];
ID(indExcludeLetter) = [];
Location(indExcludeLetter) = [];
Genotype(indExcludeLetter,:) = [];
CoorGerman = Coor(indGerman);
AncestralGerman = Ancestral(indGerman);
DerivedGerman = Derived(indGerman);

DateGermanFirstOccurance = [];
TimeGermanDAF = [];
LocationGermanFirstOccurance = [];
DAFIndices = [];
indDAFGerman = [];
DaysGerman = [];
indAnc = [];
IDGerman = [];
for i = 1:length(CoorGerman)
    indDA = find(strcmp(GenotypeGerman(:,i),DerivedGerman(i))==1);
    indAnc = [indAnc; find(strcmp(GenotypeGerman(:,i),DerivedGerman(i))==0)];
    DateGermanFirstOccurance = [DateGermanFirstOccurance; DateNumber(indDA(1:5))'];
    TimeGermanDAF = [TimeGermanDAF;MonthDay(indDA(1:5),:)'];
    LocationGermanFirstOccurance = [LocationGermanFirstOccurance; Location(indDA(1:5))'];
    DAFIndices = [DAFIndices;indDA(1:5)'];
    indDAFGerman = [indDAFGerman;{indDA}];
    DaysGerman = [DaysGerman;DateNumber(indDA(1:5))'];
    IDGerman = [IDGerman; ID(indDA(1:5))'];
end
AllGenotypeString = [];
GenotypeGermanString = [];
for i = 1:length(GenotypeGerman);
    GenotypeGermanString = [GenotypeGermanString;strcat(GenotypeGerman(i,1),GenotypeGerman(i,2),GenotypeGerman(i,3),GenotypeGerman(i,4),GenotypeGerman(i,5))];
%     strAll = {''};
%     for j = 1:27;
%         strAll = strcat(strAll,Genotype(i,j));
%     end
%     AllGenotypeString = [AllGenotypeString;strAll];
end
uniqueString = unique(GenotypeGermanString);

PathWayGerman = [{'CCCAG'};{'CCCGG'};{'TTCGG'};{'TTTGG'};{'TTTGT'}];

indAnc = find(strcmp(GenotypeGermanString,PathWayGerman{1})==1);
indStp1 = find(strcmp(GenotypeGermanString,PathWayGerman{2})==1); % the Belguin one is on a different background with additional C8782T mutation
indGerman = find(strcmp(GenotypeGermanString,PathWayGerman{3})==1); 
indGermanP1 = find(contains(GenotypeGermanString,'TTTG')==1);
indGermanP1NoLOF = find(strcmp(GenotypeGermanString,PathWayGerman{4})==1);
indGermanP1PLOF = find(strcmp(GenotypeGermanString,PathWayGerman{5})==1);
otherString = setdiff(uniqueString,PathWayGerman);
indAmerica = find(contains(Location(:,1),'America')==1);
indChina = find(contains(Location(:,1),'China')==1);
indUS = find(contains(Location(:,1),'US')==1);
indItaly = find(contains(Location(:,1),'Italy')==1);
indFrance = find(contains(Location(:,1),'France')==1);
indBelgium = find(contains(Location(:,1),'Belgium')==1);
indGermany = find(contains(Location(:,1),'Germany')==1);
indUK = find(contains(Location(:,1),'United Kingdom')==1);
indIceland = find(contains(Location(:,1),'Iceland')==1);
ResultGermany = [];
ResultFrance = [];
ResultBel = [];
ResultUK = [];
ResultIce = [];
ResultUS = [];
for d = 0:max(DateNumber)
    l1 = length(intersect(intersect(find(DateNumber <=d),indAnc),indGermany));
    %l2 = length(intersect(intersect(find(DateNumber <=d),indStp1),indEurope));
    l3 = length(intersect(intersect(find(DateNumber <=d),indGerman),indGermany));
    lN = length(intersect(intersect(find(DateNumber <=d),indGermanP1),indGermany));
    l4 = length(intersect(intersect(find(DateNumber <=d),indGermanP1NoLOF),indGermany));
    l5 = length(intersect(intersect(find(DateNumber <=d),indGermanP1PLOF),indGermany));
    ResultGermany = [ResultGermany;l1,l3,lN,l4,l5];
    l1 = length(intersect(intersect(find(DateNumber <=d),indAnc),indFrance));
    l3 = length(intersect(intersect(find(DateNumber <=d),indGerman),indFrance));
    lN = length(intersect(intersect(find(DateNumber <=d),indGermanP1),indFrance));
    l4 = length(intersect(intersect(find(DateNumber <=d),indGermanP1NoLOF),indFrance));
    l5 = length(intersect(intersect(find(DateNumber <=d),indGermanP1PLOF),indFrance));
    ResultFrance = [ResultFrance;l1,l3,lN,l4,l5];
    l1 = length(intersect(intersect(find(DateNumber <=d),indAnc),indBelgium));
    l3 = length(intersect(intersect(find(DateNumber <=d),indGerman),indBelgium));
    lN = length(intersect(intersect(find(DateNumber <=d),indGermanP1),indBelgium));
    l4 = length(intersect(intersect(find(DateNumber <=d),indGermanP1NoLOF),indBelgium));
    l5 = length(intersect(intersect(find(DateNumber <=d),indGermanP1PLOF),indBelgium));
    ResultBel = [ResultBel;l1,l3,lN,l4,l5];
    l1 = length(intersect(intersect(find(DateNumber <=d),indAnc),indUK));
    l3 = length(intersect(intersect(find(DateNumber <=d),indGerman),indUK));
    lN = length(intersect(intersect(find(DateNumber <=d),indGermanP1),indUK));
    l4 = length(intersect(intersect(find(DateNumber <=d),indGermanP1NoLOF),indUK));
    l5 = length(intersect(intersect(find(DateNumber <=d),indGermanP1PLOF),indUK));
    ResultUK = [ResultUK;l1,l3,lN,l4,l5];
    l1 = length(intersect(intersect(find(DateNumber <=d),indAnc),indIceland));
    l3 = length(intersect(intersect(find(DateNumber <=d),indGerman),indIceland));
    lN = length(intersect(intersect(find(DateNumber <=d),indGermanP1),indIceland));
    l4 = length(intersect(intersect(find(DateNumber <=d),indGermanP1NoLOF),indIceland));
    l5 = length(intersect(intersect(find(DateNumber <=d),indGermanP1PLOF),indIceland));
    ResultIce = [ResultIce;l1,l3,lN,l4,l5];
    l1 = length(intersect(intersect(find(DateNumber <=d),indAnc),indUS));
    l3 = length(intersect(intersect(find(DateNumber <=d),indGerman),indUS));
    lN = length(intersect(intersect(find(DateNumber <=d),indGermanP1),indUS));
    l4 = length(intersect(intersect(find(DateNumber <=d),indGermanP1NoLOF),indUS));
    l5 = length(intersect(intersect(find(DateNumber <=d),indGermanP1PLOF),indUS));
    ResultUS = [ResultUS;l1,l3,lN,l4,l5];
end
subplot(3,2,1);
cmap = [0.5 0.5 0.5;
    1 0.1 .1;
    0.2 0.8 0.2;
    .8 0.8 0.2;
    0.0 .3 .9];
hold on
for i = 1:5
    if i==3;
        plot([0:max(DateNumber)],ResultGermany(:,i)/sum(ResultGermany(end,:)),'--', 'LineWidth',3, 'Color', cmap(i,:));
    else
        plot([0:max(DateNumber)],ResultGermany(:,i)/sum(ResultGermany(end,:)),'-', 'LineWidth',3, 'Color', cmap(i,:));
    end
end
legend('EP-3','EP','EP+1 and EP+1+LOF','EP+1','EP+1+LOF','Location','northwest')

legend boxoff 
title('Germany','FontWeight','normal');
grid on;
box on
set(gca,'XTick',[0:30:100]);
set(gca,'FontSize',22);
xlim([25 101])
set(gca,'XTickLabel',[{'Dec-24'},{'Jan-23'},{'Feb-22'},{'Mar-23'}]);

text(25,.52,'Early founder','FontSize',25,'FontWeight','bold');
subplot(3,2,2);
hold on
for i = 1:5
    if i ==3;
        plot([0:max(DateNumber)],ResultBel(:,i)/sum(ResultBel(end,:)),'--', 'LineWidth',3, 'Color', cmap(i,:));
    else
        plot([0:max(DateNumber)],ResultBel(:,i)/sum(ResultBel(end,:)),'-', 'LineWidth',3, 'Color', cmap(i,:));
    end
end
title('Belgium','FontWeight','normal');
box on
set(gca,'XTick',[0:30:100]);
set(gca,'FontSize',22);
xlim([55 101])
set(gca,'XTickLabel',[{'Dec-24'},{'Jan-23'},{'Feb-22'},{'Mar-23'}]);
text(55,.58,'Late founder','FontSize',25,'FontWeight', 'bold');

grid on;
subplot(3,2,3);
hold on
for i = 1:5
    if i==3
    plot([0:max(DateNumber)],ResultFrance(:,i)/sum(ResultFrance(end,:)),'--', 'LineWidth',3, 'Color', cmap(i,:));
    else
            plot([0:max(DateNumber)],ResultFrance(:,i)/sum(ResultFrance(end,:)),'-', 'LineWidth',3, 'Color', cmap(i,:));
    end
end
title('France','FontWeight','normal');
set(gca,'XTick',[0:30:100]);
set(gca,'FontSize',22);
xlim([25 101])
ylabel('Cumulative frequency over time')
set(gca,'XTickLabel',[{'Dec-24'},{'Jan-23'},{'Feb-22'},{'Mar-23'}]);
grid on;
box on
subplot(3,2,4)
hold on
for i = 1:5
    if i==3
    plot([0:max(DateNumber)],ResultUK(:,i)/sum(ResultUK(end,:)),'--', 'LineWidth',3, 'Color', cmap(i,:));
    else
         plot([0:max(DateNumber)],ResultUK(:,i)/sum(ResultUK(end,:)),'-', 'LineWidth',3, 'Color', cmap(i,:));
    end
end
title('UK','FontWeight','normal');
set(gca,'XTick',[0:30:100]);
set(gca,'FontSize',22);
xlim([55 101])
set(gca,'XTickLabel',[{'Dec-24'},{'Jan-23'},{'Feb-22'},{'Mar-23'}]);
grid on;
box on
subplot(3,2,5);
hold on
for i = 1:5
    if i ==3
    plot([0:max(DateNumber)],ResultUS(:,i)/sum(ResultUS(end,:)),'--', 'LineWidth',3, 'Color', cmap(i,:));
    else
            plot([0:max(DateNumber)],ResultUS(:,i)/sum(ResultUS(end,:)),'-', 'LineWidth',3, 'Color', cmap(i,:));
    end
end

title('US','FontWeight','normal');
set(gca,'XTick',[0:30:100]);
set(gca,'FontSize',22);
xlim([25 101])
set(gca,'XTickLabel',[{'Dec-24'},{'Jan-23'},{'Feb-22'},{'Mar-23'}]);
xlabel('Time')
grid on;
box on
subplot(3,2,6);

hold on
for i = 1:5
    if i==3
    plot([0:max(DateNumber)],ResultIce(:,i)/sum(ResultIce(end,:)),'--', 'LineWidth',3, 'Color', cmap(i,:));
    else
        plot([0:max(DateNumber)],ResultIce(:,i)/sum(ResultIce(end,:)),'-', 'LineWidth',3, 'Color', cmap(i,:));
    end
end
box on
title('Iceland','FontWeight','normal');
set(gca,'XTick',[0:30:100]);
set(gca,'FontSize',22);
xlim([55 101])
set(gca,'XTickLabel',[{'Dec-24'},{'Jan-23'},{'Feb-22'},{'Mar-23'}]);
grid on;
xlabel('Time')
set(gcf,'PaperPosition',[0 0 16 12]);
saveas(1,'Fig7b.jpg');

