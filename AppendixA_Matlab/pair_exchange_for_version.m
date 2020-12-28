vec = (1:9);
for s=2:2:length(vec)
    tmp = vec(s);
    vec(s) = vec(s-1);
    vec(s-1) = tmp;
end
disp(vec)