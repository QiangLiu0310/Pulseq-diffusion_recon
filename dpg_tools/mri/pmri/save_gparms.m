function save_gparms(fname,N2)
%
% usage: save_gparms(fname,N2)
%
nsl  = length(N2);

nss = 1;
if iscell( N2{1} )
  nss  = length( N2{1} );

  N3 = cell([1 nsl*nss ]);
  for slc=1:nsl;
    for ss=1:nss;
      N3{ (slc-1)*nss + ss } = N2{slc}{ss};
    end;
  end;
  N2 = N3; clear N3;
end;

npat = length( N2{1} );


if isfield( N2{1}(1), 'n')
  sz = size( N2{1}(1).n );
else
  sz = size( N2{1}(1).coef );
end;
NN = zeros([ sz(1) sz(2) npat nsl ]);

yks = length( find( ( N2{1}(1).pattern == '*' ) + ...
                    ( N2{1}(1).pattern == '+' ) + ...
                    ( N2{1}(1).pattern == '-' ) + ...
                    ( N2{1}(1).pattern == '<' ) + ...
                    ( N2{1}(1).pattern == '>' ) ) );
xks = sz(1)/yks/sz(2);

meta = sprintf('X=%02d;Y=%02d;C=%02d;',xks,yks,sz(2));
if (nss>1)
  meta = [meta sprintf('Nsl=%02d;SMS=%d;',nsl,nss)];
end;
if isfield( N2{1}(1), 'eta')
  meta = [meta sprintf('eta=%g;',N2{1}(1).eta)];
end;
if isfield( N2{1}(1), 'chi')
  meta = [meta sprintf('chi=%g;',N2{1}(1).chi)];
end;

for c2=1:npat,
  meta = [ meta N2{1}(c2).pattern ':' ];
  for c1=1:nsl*nss,
    
    if isfield( N2{c1}(c2), 'n')
      NN(:,:,c2,c1) = N2{c1}(c2).n;
    else
      NN(:,:,c2,c1) = N2{c1}(c2).coef;
    end;

  end;
end;


fprintf([ meta '\n']);
savend(fname,NN,'dbl',meta);

