function spread(X, label)
% Plot samples of different labels with different colors.
[d,n] = size(X);
if nargin == 1
    label = ones(n,1);
end
assert(n == length(label));

color = 'brgmcyk';
m = length(color);
c = max(label);

switch d
    case 2
        view(2);
        for i = 1:c
            idc = label==i;
            scatter(X(1,idc),X(2,idc),36,color(mod(i-1,m)+1));
        end
    case 3
        view(3);
        for i = 1:c
            idc = label==i;
            scatter3(X(1,idc),X(2,idc),X(3,idc),36,color(mod(i-1,m)+1)); hold on;
        end
        hold on; [x,y,z] = sphere; 
        surface(x,y,z,'FaceColor', 'none'); hold off;
    otherwise
        error('ERROR: only support data of 2D or 3D.');
end
axis equal
grid on
hold off