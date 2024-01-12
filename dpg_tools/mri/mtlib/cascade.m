function h = cascade(IO,mm,nn,ca)
% CASCADE show a tiled version of the image sequence 
% (code borrowed from montage.m)
% usage: cascade(IO,mm,nn,ca)
%        IO -  3D image sequence to display
% optional args:
%        mm -  number of images down
%        nn -  number of images across
%        ca -  color axis

if length(size(IO)) == 2,
  IO = reshape(IO, sqrt(size(IO,1)), sqrt(size(IO,1)), size(IO,2) );
end;

siz = [size(IO,1) size(IO,2) size(IO,3)];
if nargin < 2,
  nn = sqrt(prod(siz))/siz(2);
  mm = siz(3)/nn;
end;

if nargin==2
  if length(mm) > 2; ca = mm(3); end;
  if length(mm) > 1; nn = mm(2); end;
  mm = mm(1);
end;

if (ceil(nn)-nn) < (ceil(mm)-mm),
  nn = ceil(nn); mm = ceil(siz(3)/nn);
else
  mm = ceil(mm); nn = ceil(siz(3)/mm);
end
  
  
b = IO(1,1); % to inherit type 
b(1,1) = 0;  % from a
b = zeros( mm*siz(1), nn*siz(2) );

rows = 1:siz(1); cols = 1:siz(2);
for i=0:mm-1,
  for j=0:nn-1,
    k = j+i*nn+1;
    if k<=siz(3),
      b(rows+i*siz(1),cols+j*siz(2)) = IO(:,:,k);
    end
  end
end

if nargout == 0,
imagesc(real(b)); axis('image','off');
if exist('ca'), caxis(ca); end;
%ah = colorbar;
%pos = get(ah,'position');
%set(ah,'position', [pos(1:2)  pos(3)/2 pos(4) ]);

% set(gcf,'PaperPosition',[0 0 5.5 8]);
C = colormap;
i = siz(1)*[0:mm] + 0.5;
j = siz(2)*[0:nn] + 0.5;

%% draw a line to seperate the rows:
line( [ zeros(size(i)); size(b,2)*ones(size(i)) ], [ ones(2,1)*i ], ...
     'color', C(size(C,1),:) );

%% draw a line to seperate the columns:
line( [ ones(2,1)*j ], [ zeros(size(j)); size(b,1)*ones(size(j)) ], ...
     'color', C(size(C,1),:) );
end;

if nargout > 0,
  h = b;
end;