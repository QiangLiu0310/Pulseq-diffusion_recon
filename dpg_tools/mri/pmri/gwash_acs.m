%% GRAPPA-wash the ACS data...
if ~exist('perslc'); perslc = 1:hdr.NSlc; end;
if ~exist('verbose','var'), verbose=0; end;
sz = size(A);

clear tmpA tmpB tmpC A3 B3 C3

fprintf('%d',mod(1:length(perslc),10)); fprintf('\n');
for slc = perslc; % 1:hdr.NSlc;
  fprintf('.');
  [~,~,Ng{slc}]=recongrappa_multik(sz([2 1 4]),permute(C(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R);

  for cnt=1:hdr.R;
    Fa(:,:,:,cnt) = recongrappa_multik(sz([2 1 4]),permute(A(:,cnt:hdr.R:end,slc,:),[2 1 4 3]),cnt:hdr.R:sz(2),'kernel','2x5','N',Ng{slc});
    Fb(:,:,:,cnt) = recongrappa_multik(sz([2 1 4]),permute(B(:,cnt:hdr.R:end,slc,:),[2 1 4 3]),cnt:hdr.R:sz(2),'kernel','2x5','N',Ng{slc});
  end;

  tmpA(:,:,slc,:) = mean(Fa,4);
  tmpB(:,:,slc,:) = mean(Fb,4);
  tmpC(:,:,slc,:) = ifi( phzshift( fif(tmpA(:,:,slc,:)), fif(tmpB(:,:,slc,:)), {'nofft'}) );
  [~,~,Ng2{slc}]=recongrappa_multik(sz([2 1 4]),tmpC(:,:,slc,:),[],'kernel','2x5','dks',hdr.R*[1 2]);

  for cnt=1:2*hdr.R;
    Fa(:,:,:,cnt) = recongrappa_multik(sz([2 1 4]),permute(A(:,cnt:2*hdr.R:end,slc,:),[2 1 4 3]),cnt:2*hdr.R:sz(2),'kernel','2x5','N',Ng2{slc});
    Fb(:,:,:,cnt) = recongrappa_multik(sz([2 1 4]),permute(B(:,cnt:2*hdr.R:end,slc,:),[2 1 4 3]),cnt:2*hdr.R:sz(2),'kernel','2x5','N',Ng2{slc});
  end;

  A3(:,:,slc,:) = mean(Fa,4);
  B3(:,:,slc,:) = mean(Fb,4);
  C3(:,:,slc,:) = ifi( phzshift( fif(A3(:,:,slc,:)), fif(B3(:,:,slc,:)), {'nofft'}) );

  if verbose; keyboard; end;
end;
fprintf('\n');

A3 = permute(A3,[2 1 3 4]);
B3 = permute(B3,[2 1 3 4]);
B3 = permute(C3,[2 1 3 4]);

A4 = zeros(size(A3));
B4 = zeros(size(B3));
C4 = zeros(size(C3));
for slc=perslc; fprintf('.'); for coil=1:size(A3,4); for cnt=1:2*hdr.R; 
 A4(:,cnt:2*hdr.R:end,slc,coil) = iffts(fft(iffts(phzshift( ffts(ifft(ffts(A3(:,cnt:2*hdr.R:end,slc,coil).',2),[],2),2), ffts(ifft(ffts(A(:,cnt:2*hdr.R:end,slc,coil).',2),[],2),2),{'nofft','nocombo'}),2),[],2),2).'; 
 B4(:,cnt:2*hdr.R:end,slc,coil) = iffts(fft(iffts(phzshift( ffts(ifft(ffts(B3(:,cnt:2*hdr.R:end,slc,coil).',2),[],2),2), ffts(ifft(ffts(B(:,cnt:2*hdr.R:end,slc,coil).',2),[],2),2),{'nofft','nocombo'}),2),[],2),2).'; 
end; end; end; fprintf('\n');

for slc=perslc;
  C4(:,:,slc,:) = ifi( permute(phzshift( fif(permute(A4(:,:,slc,:),[2 1 3 4])), ...
                                         fif(permute(B4(:,:,slc,:),[2 1 3 4])), {'nofft'} ),[2 1 3 4]) );
end;

A3(:,:,hdr.slcindx(perslc),:) = A3(:,:,perslc,:);
B3(:,:,hdr.slcindx(perslc),:) = B3(:,:,perslc,:);
C3(:,:,hdr.slcindx(perslc),:) = C3(:,:,perslc,:);

A4(:,:,hdr.slcindx(perslc),:) = A4(:,:,perslc,:);
B4(:,:,hdr.slcindx(perslc),:) = B4(:,:,perslc,:);
C4(:,:,hdr.slcindx(perslc),:) = C4(:,:,perslc,:);
