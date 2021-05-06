K = 16000;
nbins = 32000;
a = linspace(3.5,4,K+1);
diagram = 255*ones(nbins,K+1,'uint8');

N = 10000;
M = 10000;

for i = 1:K+1
    if a(i) ~= 4
        y = 1/2;
    else
        y = 1/3;
    end
    for j = 1:N
        y = a(i)*y*(1-y);
    end
    for j = 1:M
        y = a(i)*y*(1-y);
        idx = floor(y*nbins)+1;
        diagram(idx,i) = 0;
    end
end

diagram = flipud(diagram);
img = imresize(diagram,1/8,'bicubic');
imshow(img);
imwrite(img,'../../figures/capitolo7/bifurcation-diagram-zoomed.png');