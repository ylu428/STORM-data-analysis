function comp = GFPfilter(data, MaxInt) % delete all of the data with intensity above "MaxInt"
    comp = data;
    for i = 1:length(data)
        if data(i,6) >= MaxInt
            comp(i,1:5) = NaN;
        end
    end

end