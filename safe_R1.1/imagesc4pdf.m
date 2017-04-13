%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Plots the pattern plot fir the sensitivity indices      %%%%%%

% to first line of zeros for O2 SWI flux Index:
%SIndex400m_zeros = [zeros(1, size(SIndex400m,2),1); SIndex400m; zeros(1, size(SIndex400m,2),1)];
% combine 400m and 400m and use as input for this function
%SIndex_400_4000m = [SIndex400m_zeros; SIndex4000m];

function h = imagesc4pdf(C)
labelparams = { 'k1','k2ord','f1','KNH4','gamma NH4','gamma H2S','K I','K II','ks' 'km' 'ka' } ; % input names
Titles = {'O_2', 'NO_3', 'SO_4', 'NH_4', 'H_2S', 'PO_4', 'DIC', 'ALK', 'zeros'};


[ny nx] = size(C);

px = bsxfun(@plus, [-0.5; 0.5; 0.5; -0.5], reshape(1:nx, [1 1 nx]));
py = bsxfun(@plus, [-0.5; -0.5; 0.5; 0.5], 1:ny);

n = numel(C);
px = reshape(repmat(px, [1 ny 1]), 4, n);
py = reshape(repmat(py, [1 1 nx]), 4, n);

h = patch(px, py, reshape(C,1,n), 'linestyle', '-');

%set(0,'defaulttextinterpreter','latex')

xlim([.5 nx+.5]);
ylim([.5 ny+.5]);
set(gca, 'ydir', 'reverse');
set(gca,'YTick',[1:1:ny]);
set(gca, 'yticklabel',Titles);
set(gca, 'xticklabel',labelparams);
colormap autumn
colormap(flipud(colormap))
colorbar
print('-depsc2', ['0_KSIndex_ALL_OUTPUT_13042017.eps']);