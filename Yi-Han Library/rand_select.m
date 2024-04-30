function output = rand_select(input, target, current)
    arguments
        input
        target double =1;
        current double =1;
            end
    rng(1)
    y = rand(size(input,1),1);
    rati = target/current;
    output = input(y<rati,:);
end

