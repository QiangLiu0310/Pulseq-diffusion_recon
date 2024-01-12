function N = load_gparms( fname )

[g,m]=readnd( fname );
str = sprintf('%c',m);
fprintf('grappa desc: %s\n',str);
indx = [ findstr(str,';') findstr(str,':') findstr(str,'#') ]-1;

SMS = 1;
% evaluate all of the parameter declarations in the string:
t = regexp(str,'(\w*)=([\d\.e\+\-]*);','tokens');
for cnt=1:length(t);
  % fprintf([ t{cnt}{1} ' = ' t{cnt}{2} ';\n']); 
  eval([ t{cnt}{1} ' = ' t{cnt}{2} ';']); 
end;

N = {};
cnt0 = 1; pcnt = 1;
for cnt=1:length(indx)
  tmp = str(cnt0:indx(cnt));
  if ( length([ findstr( tmp, 'x' ) findstr(tmp,'>') findstr(tmp,'<') ]) == 1 ),
    dk = max(diff( findstr(tmp,'*') ));
    if ( length(findstr(tmp,'*')) == 1 )
      dk = 1;
    end;
    pattern = tmp;
    for Scnt=1:size(g,4)
     N{Scnt}(pcnt,1).dk = dk;
     N{Scnt}(pcnt,1).pattern = tmp;
     N{Scnt}(pcnt,1).coef = g(:,:,pcnt,Scnt);
     if (dk == 1)
       if ( str(indx(cnt)+1) == '#' )
         N{Scnt}(pcnt,1).dir = 'parallel';
       else,
         N{Scnt}(pcnt,1).dir = 'normal';
       end;
     end;
    end;

    pcnt = pcnt + 1;
  end;
  cnt0 = indx(cnt)+2;
end;

if (SMS>1)
  for slc=1:Nsl;
    for cnt=1:SMS;
      N2{slc}{cnt} = N{ (slc-1)*SMS + cnt };
    end;
  end;
  N = N2; clear N2;

  if exist('eta','var'), N{1}{1}(1).eta = eta; end;
  if exist('chi','var'), N{1}{1}(1).chi = chi; end;
else
  if exist('eta','var'), N{1}(1).eta = eta; end;
  if exist('chi','var'), N{1}(1).chi = chi; end;  
end;

return;
