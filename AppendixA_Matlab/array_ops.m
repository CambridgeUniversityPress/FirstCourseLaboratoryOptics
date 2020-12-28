% This code is intended to be added line by line to the console/command window
A = [ [2,0,1]; [0,1,3]; [-1,1,0] ]
B = [ [1,1,2]; [1,2,-2]; [1,4,1] ]
A(2,1)       % single element
A(2,:)       % entire row
A(:,2)       % entire column
A.'          % transpose
A*B          % matrix mult.
B*A			
A.*B         % element-wise mult.