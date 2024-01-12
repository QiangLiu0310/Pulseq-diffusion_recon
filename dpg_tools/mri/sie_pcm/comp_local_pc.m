function phz = comp_local_pc( S1, S2, verbose )
%%
%% an implementation of the "local phase correction" that is the new
%% Siemens standard.  Based on the patent by Feiweier and JonP's code.
%%

% S1 = ffts(ifft(ffts( kin{1}(64,:,1).' ,1),[],1),1);
% S2 = ffts(ifft(ffts( kin{2}(64,:,1).' ,1),[],1),1);

if nargin<3, verbose=0; end;

[m c]= size(S1);

indx = 1:(m-1);

for cnt=1:c,                            % for each coil...
  
  % compute the cross-correlation between the two navigator signals
  C0 = S2(indx+0,cnt).*conj(S1(indx+0,cnt));
  C1 = S2(indx+1,cnt).*conj(S1(indx+1,cnt));

  % 'average' the phase difference between them.  This yeilds the linear
  % phase term:
  dphi1 = angle( sum( C1 .* conj(C0), 1 ) );
  
  % apply the linear phase term to the phase difference signals,
  indx2 = [1:m]-m/2;
  cr = S1(:,cnt) .* conj(S2(:,cnt)) .* exp(j*indx2(:)*dphi1);

  % use the phase average here to find the constant phase term.
  dphi0 = angle(sum( cr ,1));
  
  % return the data, in a format ready for phzapply.m
  phz(:,cnt) = [ dphi1;  ...
                 dphi0 ]; 

end;

if (verbose)
  pltcmplx( abs(S1(:,:)).*exp(j*angle(S1(:,:)./S2(:,:))), abs(S2(:,:)).*exp(j*(ones(size(indx2(:)))*phz(2,:) - (indx2(:))*phz(1,:)) ) );
end;

sz = size(S1);
if length(sz)>2
  phz = reshape(phz,[2 sz(2:end)]);
end;

