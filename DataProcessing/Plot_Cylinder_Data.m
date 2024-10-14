data = readtable("Levitation Table of Levitation - Sheet1.csv");

material = data.Material;
length_mm = data.Length_mm_;
diameter_um = data.Diameter__m_;
maxpressure_torr = data.MaxPressure;
minpressure_torr = data.MinPressure;
passfail = data.P_F;

Nsamples = length(material);

matlist = {};
legendentries = [];
legendlabels = {};
colorlist = colororder();
markerlist = {"-o", "-square", "-diamond", "-v", "-<", "->", "-*"};
addtolegend = false;

figure();
hold on;
for n = 1:Nsamples
    % check if we've already plotted a material of this type 
    midx = find(strcmp(matlist, material{n}));

    if isempty(midx) 
        if isempty(material{n})
            continue  % skip entries with no material listed
        end
        % if we haven't had a material like this yet
        matlist{end+1} = material{n}; % add the material to the list
        midx = length(matlist); % update materialindex
        addtolegend = true; % flag to add this line to the legend
    end
    % plot the line 
    p = plot([diameter_um(n), diameter_um(n)], ...
        [minpressure_torr(n), maxpressure_torr(n)], ...
        markerlist{midx}, "Color", colorlist(midx, :));

    % if this is the first material of this type, add to legend 
    if addtolegend
        legendentries(midx) = p;
        legendlabels{midx} = matlist{midx};
    end
end
legend(legendentries, legendlabels);
xlabel("Diameter [µm]");
ylabel("Levitation Pressure Range [Torr]")

% this is just to widen the plot range a bit 
scmat = [1.1, -0.1; -0.1, 1.1]; 
xl = xlim();
yl = ylim();
xlim(xl*scmat);
ylim(yl*scmat);
