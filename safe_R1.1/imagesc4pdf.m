function h = imagesc4pdf(C)

[ny nx] = size(C);

px = bsxfun(@plus, [-0.5; 0.5; 0.5; -0.5], reshape(1:nx, [1 1 nx]));
py = bsxfun(@plus, [-0.5; -0.5; 0.5; 0.5], 1:ny);

n = numel(C);
px = reshape(repmat(px, [1 ny 1]), 4, n);
py = reshape(repmat(py, [1 1 nx]), 4, n);

h = patch(px, py, reshape(C,1,n), 'linestyle', '-');

xlim([.5 nx+.5]);
ylim([.5 ny+.5]);
set(gca, 'ydir', 'reverse');
colormap copper
colormap(flipud(colormap))
colorbar
print('-depsc2', ['RESULTS/SIndex_500m.eps']);