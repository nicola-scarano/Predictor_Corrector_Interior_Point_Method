%% Function that generates the matrix A 
%This function creates a matrix A having: the first row with ones in the
%even positions and zeros in the odd positions; the second row with ones
%in the odd positions and zeros in the even positions; the third row with
%ones in the first half and minus 1 in the second half; 
%It returns also the vector b as a all one column vector with size (3,1)
%and the vector c

function [A, b, c] = Abc_setup(n, a)
     
% generation of A
    A = ones(3,n);

    A(1,1:2:n) = 0;

    A(2,2:2:n) = 0;

    A(3,n/2+1:1:n) = -1;

%   generation of b
    b = [1;1;0];

%   generation of c
%generation of c
    c = zeros(n,1);
    for i = 1:n
        l = log(i);
        if(mod(i,2) == 0)
            c(i) = l*a;
        else
            c(i) = l;
        end
    end
end
    
    
    
