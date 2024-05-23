function MolCell2021_plotfig1(~)
with_boot_figs = 0;
h1 = figure;

DhistMinSteps = 4;
Dfit_inital_guess = [0.2 0.11 0.8 0.31];

file_path = '..\Mendeley_Data';
protein_names = {'MukB', 'RNAP', 'GyrA', 'ParC', 'UvrA', 'MutS', 'LacI', 'Pol1', 'LigA', 'HNS', 'HU'};

nprots = length(protein_names);
number_of_cells_all = zeros(nprots,1);
prot_frac = zeros(nprots,1);
prot_frac_err = zeros(nprots,2);
prot_dmob = zeros(nprots,1);
prot_dmob_err = zeros(nprots,2);
number_of_tracks = zeros(nprots,1);
prot_dmob_single = zeros(nprots,1);
for prot = 1:nprots
    %%

    rangeD = 0:0.02:2.5;
    binSpacing = rangeD(2)-rangeD(1);
    x = rangeD+binSpacing;
    
    folder_names = dir([file_path '\' protein_names{prot} '\*.mat']);
    
    D_all_cells = zeros(1e6,1) ;
    cell_index = zeros(1e6,1) ;
    cell_index_temp = 1;
    kk = 1;
    for jj = 1:length(folder_names)
        load([file_path '\' protein_names{prot} '\' folder_names(jj).name]);
        for ll = 1:length(PALMcell)
            
            if ~isempty(PALMcell(ll).tracks)
                D = tracks_to_D(PALMcell(ll).tracks);
                if ~isempty(D)
                    D_all_cells(kk:kk+length(D)-1,1) = D;
                    cell_index(kk:kk+length(D)-1,1) = cell_index_temp*ones(length(D),1);
                    cell_index_temp = cell_index_temp+1;
                    kk = kk+length(D);
                end
            end
        end
        
    end
    D_all_cells(kk:end) = [];
    cell_index(kk:end) = [];
    number_of_cells = max(cell_index);
    number_of_tracks(prot) = kk -1;
    
    figure(h1)
    subplot(4,3,prot)
    histDCount = histc(D_all_cells,rangeD + binSpacing/2);
    histDCount = histDCount./(sum(histDCount)*binSpacing);
    hold on
    bar1 = bar(rangeD+binSpacing,histDCount,'BarWidth',1,'FaceColor',[.5 .5 .5],'LineStyle','none');
    xlim([min(rangeD), max(rangeD)]);
    text(1, 1, [num2str(number_of_cells) ' cells ' num2str(number_of_tracks(prot)) ' mols'])
    
    [f, mob] = Dhist_fitting(histDCount,x',DhistMinSteps, Dfit_inital_guess );
    
    hold on;
    x2 = linspace(0,max(x),1000);
    cFit1 = one_gamma_fit([f/100,Dfit_inital_guess(2)],x2,DhistMinSteps);
    plot(x2, cFit1, 'r-','LineWidth',1)
    
    hold on;
    cFit2 = one_gamma_fit([(1-f/100),mob],x2,DhistMinSteps);
    plot(x2, cFit2, 'b-','LineWidth',1)
    
    title({protein_names{prot}, [num2str(f,3), '% at D = ', num2str(Dfit_inital_guess(2),2)],[num2str(100-f,3), '% at D = ', num2str(mob,2)] })
    
    prot_dmob_single(prot) = mob;
    
    
    
    %% bootstrap cells
    nboots = 100;
    frac1 = zeros(1,nboots);
    D_mob = zeros(1,nboots);
    %prot
    %number_of_cells
    if with_boot_figs == 1
        h2 = figure;
    end
    for iboot = 1:nboots
        
        y = datasample(1:number_of_cells,number_of_cells);
        D_to_use = zeros(1e6, 1);
        jj = 1;
        for ii = 1:length(y)
            tempD =    D_all_cells(cell_index == y(ii));
            D_to_use(jj:jj+length(tempD)-1,1) = tempD;
            jj = jj+length(tempD);
        end
        D_to_use(jj:end) = [];
        histDCount = histc(D_to_use,rangeD + binSpacing/2);
        histDCount = histDCount./(sum(histDCount)*binSpacing);
        
        if with_boot_figs == 1
            figure(h2)
            clf
            bar1 = bar(rangeD+binSpacing,histDCount,'BarWidth',1,'FaceColor',[.5 .5 .5],'LineStyle','none');
            xlim([min(rangeD), max(rangeD)]);
                       
            [f, mob] = Dhist_fitting(histDCount,x',DhistMinSteps, Dfit_inital_guess );
            frac1(iboot) = f;
            D_mob(iboot) = mob;
            
            hold on;
             x2 = linspace(0,max(x),1000);
             cFit1 = one_gamma_fit([frac1(iboot)/100,Dfit_inital_guess(2)],x2,DhistMinSteps);
             plot(x2, cFit1, 'r-','LineWidth',1)
    
                hold on;
             cFit2 = one_gamma_fit([(1-frac1(iboot)/100),D_mob(iboot)],x2,DhistMinSteps);
             plot(x2, cFit2, 'b-','LineWidth',1)
    
            title({protein_names{prot}, [num2str(frac1(iboot),3), '% at D = ', num2str(Dfit_inital_guess(2),2)],[num2str(100-frac1(iboot),3), '% at D = ', num2str(D_mob(iboot),2)] })
    
            pause(1)
            
        else
            [f, mob] = Dhist_fitting(histDCount,x',DhistMinSteps, Dfit_inital_guess );
             frac1(iboot) = f;
            D_mob(iboot) = mob;
           
        end
        
    end
    
    
    number_of_cells_all(prot) =number_of_cells;
    prot_frac(prot) = mean(frac1);
    prot_frac_err(prot,:) = prctile(frac1,[2.5 97.5]);
    prot_dmob(prot) = mean(D_mob);
    prot_dmob_err(prot,:) = prctile(D_mob,[2.5 97.5]);
    
end
formatSpec = '%.2f [%.2f - %.2f]\n';
disp('Dmob:')
for mm= 1:length(prot_frac)
        fprintf(formatSpec,prot_dmob(mm),prot_dmob_err(mm,1), prot_dmob_err(mm,2));
        Dmob_table1{mm,1}= sprintf(formatSpec,prot_dmob(mm),prot_dmob_err(mm,1), prot_dmob_err(mm,2));
end
formatSpec = '%.0f [%.0f - %.0f]\n';
disp('Frac bound:')
for mm= 1:length(prot_frac)
        fprintf(formatSpec,prot_frac(mm),prot_frac_err(mm,1), prot_frac_err(mm,2))
        frac_bound_table1{mm,1}= sprintf(formatSpec,prot_frac(mm),prot_frac_err(mm,1), prot_frac_err(mm,2));
end
figure
set(gcf,'color','w')

table1 = table(protein_names', Dmob_table1, frac_bound_table1, 'VariableNames',{'Protein','Dmob','fracbound'});
filename = 'Table1.xlsx';
writetable(table1,filename)
error_l =[prot_frac-prot_frac_err(:,1)];
error_u =[prot_frac-prot_frac_err(:,2)];
b1= bar(prot_frac);
hold on
errorbar(b1.XData',b1.YData',error_l, error_u,'k', 'LineStyle','none')
xticks(1:11)
xticklabels(protein_names)
xtickangle(45)
ylabel('% long lived binding')
figure
set(gcf,'color','w')

error_l =[prot_dmob-prot_dmob_err(:,1)];
error_u =[prot_dmob-prot_dmob_err(:,2)];
b1= bar(prot_dmob);
hold on
errorbar(b1.XData',b1.YData',error_l, error_u,'k', 'LineStyle','none')

xticklabels(protein_names)
xtickangle(45)
ylabel('mobility ')

%numbers_fortable_1

%% Fig S1
copy_number = [100 4000 80 500 800 100 80 100 40 20000 50000 80000]

figure
hold on
for ii = 1:length(prot_dmob_single)
scatter(copy_number(ii), prot_dmob_single(ii) )
end
set(gca, 'XScale', 'log')
legend(protein_names)
xlabel('literature copy number')
ylabel('measured D_mobile')
axis([10 1e6 0 3])
end
%%
function D = tracks_to_D(tracks)
dT = 0.01548;
DhistMinSteps = 4;
pixel = 0.096;
single_frame_steps= find((tracks(2:end,3) - tracks(1:end-1,3))==1);
%find elements in tracks from the same molecules as consective element
same_mol_steps= find((tracks(2:end,4) - tracks(1:end-1,4))==0);
%intersect of both gives the single frame steps only from the same molecule
single_frame_same_mol = intersect( same_mol_steps, single_frame_steps);

%Calculate teh squared dispacement form these elements in tracks
single_frame_sq_displacements = sum((tracks(single_frame_same_mol,1:2) - tracks(single_frame_same_mol+1,1:2)).^2,2);

%molecule numbers from only the desired track elements
single_frame_tracks = tracks(single_frame_same_mol,4);
% histogram of mol numbers gives the step number of each desired track
n = histc(single_frame_tracks,1:max(single_frame_tracks));
%find tracks with  enough steps
long_enough_tracks = find(n>=DhistMinSteps);

%pre allocate variables
MSD_all = zeros(length(long_enough_tracks),1);

for ii = 1:length(long_enough_tracks)
    
    %find the elements of the track
    xx = find(single_frame_tracks == long_enough_tracks(ii));
    
    %calulate the mean of the squared displacements and store values in
    %array. Truncated tracks so only take first DhistMinSteps steps
    MSD_all(ii,1) = mean(single_frame_sq_displacements(xx(1:DhistMinSteps)));
    
end
MSD = MSD_all * pixel^2;
D = MSD/(4*dT);

end

%%
function [frac1, D_mob] = Dhist_fitting(histDCount,x,n, initGuess)
curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions, 'Display', 'off', 'MaxFunEvals', 100000, 'MaxIter', 100000);
x2 = linspace(0,max(x),1000);

initGuess2 = [initGuess(1) initGuess(4)];
constrained_fit_vals = initGuess(2);
lb = [0 0];
ub = [inf inf];

c  = lsqcurvefit(@(parameters,x) ...
    two_gamma_1constrained_fit(...
    parameters,x,n,constrained_fit_vals),...
    initGuess2, x, histDCount, lb, ub, curvefitoptions);

D_mob = c(2);
frac1 = c(1)*100;
frac2 = 100 - frac1;

end

function F = two_gamma_1constrained_fit( parameters, x, n, b)

A(1) = parameters(1);
D(1) = parameters(2);

F= A(1)*(n/b)^n*x.^(n-1).*exp(-n*x/b)/factorial(n-1) + ...
    (1 - A(1))*(n/D(1))^n*x.^(n-1).*exp(-n*x/D(1))/factorial(n-1);


end

function F = one_gamma_fit( parameters, x, n)

A(1) = parameters(1);
D(1) = parameters(2);

F= A(1)*(n/D(1))^n*x.^(n-1).*exp(-n*x/D(1))/factorial(n-1);

end