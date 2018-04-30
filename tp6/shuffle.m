function v = shuffle(a)
    v = a(randperm(length(a)),:);
end
