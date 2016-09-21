function h = imagesc4pdf(C)
Titles = {'O_2', 'NO_3', 'SO_4', 'NH_4', 'H_2S', 'PO_4'};

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
%set(gca, 'yticklabel',{'$NO_3$' 'SO_4' 'NH_4' 'H_2S' 'null' 'O_2' 'NO_3' 'SO_4' 'NH_4' 'H_2S'},'Interpreter','LaTex');
set(gca, 'xticklabel',[]);
colormap autumn
colormap(flipud(colormap))
colorbar
print('-depsc2', ['RESULTS_PAWN/KSIndex_ALL_autumn.eps']);