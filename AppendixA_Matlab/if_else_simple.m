% a = int8(10);  % gives integer arithmetic
a = 10;
b = 3;

if a/b==3
    disp(' ')
    disp('Matlab is doing integer arithmetic.')
else
    disp(' ')
    disp('Matlab is doing floating point arithmetic.')
end