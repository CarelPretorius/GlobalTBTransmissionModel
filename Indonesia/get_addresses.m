function [i, s, d, lim] = get_addresses(groups, i, s, d, lim)

% Initiate any sets not so far covered by s
if ~isempty(s)
    fnames = fieldnames(s);
else
    fnames = {};
end
for ig = 1:length(groups)
    gp = groups{ig};
    for ig2 = 1:length(gp)
        if ~ismember(gp{ig2},fnames)
            s.(gp{ig2}) = [];
        end
    end
end

if length(groups) == 1
    gp1 = groups{1};
    for ig1 = 1:length(gp1)
        lim = lim+1;
        i.(gp1{ig1}) = lim;
        s.(gp1{ig1}) = [s.(gp1{ig1}), lim];
        d{lim} = [gp1{ig1}];
    end
end

if length(groups) == 2
    gp1 = groups{1}; gp2 = groups{2};
    
    for ig1 = 1:length(gp1)
        for ig2 = 1:length(gp2)
            lim = lim+1;
            i.(gp1{ig1}).(gp2{ig2}) = lim;
            s.(gp1{ig1}) = [s.(gp1{ig1}), lim];
            s.(gp2{ig2}) = [s.(gp2{ig2}), lim];
            d{lim} = [gp1{ig1}, ' ', gp2{ig2}];
        end
    end
end

if length(groups) == 3
    gp1 = groups{1}; gp2 = groups{2}; gp3 = groups{3};
    for ig1 = 1:length(gp1)
        for ig2 = 1:length(gp2)
            for ig3 = 1:length(gp3)
                lim = lim+1;
                i.(gp1{ig1}).(gp2{ig2}).(gp3{ig3}) = lim;
                s.(gp1{ig1}) = [s.(gp1{ig1}), lim];
                s.(gp2{ig2}) = [s.(gp2{ig2}), lim];
                s.(gp3{ig3}) = [s.(gp3{ig3}), lim];
                d{lim} = [gp1{ig1},  ' ', gp2{ig2},  ' ', gp3{ig3}];
            end
        end
    end
end

if length(groups) == 4
    gp1 = groups{1}; gp2 = groups{2}; gp3 = groups{3}; gp4 = groups{4};
    for ig1 = 1:length(gp1)
        for ig2 = 1:length(gp2)
            for ig3 = 1:length(gp3)
                for ig4 = 1:length(gp4)
                    lim = lim+1;
                    i.(gp1{ig1}).(gp2{ig2}).(gp3{ig3}).(gp4{ig4}) = lim;
                    s.(gp1{ig1}) = [s.(gp1{ig1}), lim];
                    s.(gp2{ig2}) = [s.(gp2{ig2}), lim];
                    s.(gp3{ig3}) = [s.(gp3{ig3}), lim];
                    s.(gp4{ig4}) = [s.(gp4{ig4}), lim];
                    d{lim} = [gp1{ig1},  ' ', gp2{ig2},  ' ', gp3{ig3},  ' ', gp4{ig4}];
                end
            end
        end
    end
end


if length(groups) == 5
    gp1 = groups{1}; gp2 = groups{2}; gp3 = groups{3}; gp4 = groups{4}; gp5 = groups{5};
    for ig1 = 1:length(gp1)
        for ig2 = 1:length(gp2)
            for ig3 = 1:length(gp3)
                for ig4 = 1:length(gp4)
                    for ig5 = 1:length(gp5)
                    lim = lim+1;
                    i.(gp1{ig1}).(gp2{ig2}).(gp3{ig3}).(gp4{ig4}).(gp5{ig5}) = lim;
                    s.(gp1{ig1}) = [s.(gp1{ig1}), lim];
                    s.(gp2{ig2}) = [s.(gp2{ig2}), lim];
                    s.(gp3{ig3}) = [s.(gp3{ig3}), lim];
                    s.(gp4{ig4}) = [s.(gp4{ig4}), lim];
                    s.(gp5{ig5}) = [s.(gp5{ig5}), lim];
                    d{lim} = [gp1{ig1},  ' ', gp2{ig2},  ' ', gp3{ig3},  ' ', gp4{ig4},  ' ', gp5{ig5}];
                    end
                end
            end
        end
    end
end
i.nstates = lim;
