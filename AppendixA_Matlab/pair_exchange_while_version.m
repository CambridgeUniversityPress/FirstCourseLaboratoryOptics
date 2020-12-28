vec = (1:9);
s = 2;
while s <= length(vec)
    tmp = vec(s);
    vec(s) = vec(s-1);
    vec(s-1) = tmp;
    s = s + 2;
end
disp(vec)