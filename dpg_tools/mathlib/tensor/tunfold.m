function result = tunfold(A,d)
%
% TUNFOLD - tensor unfold - reduces the dimension of a tensor by "unfolding" along one direction
%          (performs unfolding as described by deLathauwer in SIMAX 21(4):1253.)
%
% Example: given a 3D tensor, A, 
%
% unfold(A,2) =>        I1
%                     +----+ +----+ +----+
%                  I2 |    | |    | |    |
%                     |    | |    | |    |
%                     +----+ +----+ +----+
%                        \ ____|____ /
%                              I3
%
% for a 3D tensor, 
%
% unfold(A,1) =  [ A(:,1,:)  A(:,2,:)  ... A(:,n2,:)  ],
% unfold(A,2) = ( A(:,:,1)' A(:,:,2)' ... A(:,:,n3)' ),
% unfold(A,3) = ( A(1,:,:)' A(2,:,:)' ... A(n1,:,:)' )
%
% To fold the matrix back (i.e. undo the unfold), one *must* know the
% size of the final matrix and use the same permutation as the unfolding.  
% An example:
%    sz = size(A);
%    A == permute( reshape( unfold(A,d) , sz([ d:-1:1 N:-1:(d+1) ]) ), ...
%                  [  d:-1:1 N:-1:(d+1) ] );
%

% Copyright 2002, W Scott Hoge (shoge at ieee dot org) 

sz = size(A);
N = length(sz);
if (d > length(sz)), result = A; return; end;

result = reshape( permute( A, [ d:-1:1 N:-1:(d+1) ]), ...
                  [ sz(d) prod(sz([1:(d-1) (d+1):N])) ] );

% unfold(A,1) = reshape( permute( A, [ 1 3 2 ]), [sz(1) prod(sz(3:2)) ] )
% unfold(A,2) = reshape( permute( A, [ 2 1 3 ]), [sz(2) prod(sz([1 3])) ] )
% unfold(A,3) = reshape( permute( A, [ 3 2 1 ]), [sz(3) prod(sz([1 2])) ] )

