function h = imagesc4pdf(C)
labelparams = {'a', 'nu', 'w0', 'Db0', 'OC0', 'beta', 'por0'} ; % input names
Titles = {'F1', 'FF10', 'F2', 'FF100', 'F3', 'FF1000', 'F4', 'FF10000', 'dF1', 'dFF10', 'dF2', 'dFF100', 'dF3', 'dFF1000', 'dF4', 'dFF10000'};

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
print('-depsc2', ['0_KSIndex_ALL_autumn.eps']);