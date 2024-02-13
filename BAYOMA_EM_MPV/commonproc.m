function out = commonproc(in,fun)
% process multi-band using single band procedures
% This function is used by MODEIDFFT00, MODEIDFFT02, MODEIDFFT99
% to perform modeid for multi-bands using single band algorithms

% 141011 - added checking nseg

%==========================================================================
% check mandatory fields from IN
if ~isfield(in,'tdata')
  error('tdata is required field.');
end
if ~isfield(in,'fs')
  error('fs is required field');
end
if ~isfield(in,'f0')
  error('f0 is required field.');
end
if ~isfield(in,'f1f2')
  error('f1f2 is required field.');
end

%==========================================================================
% check legitimacy

% ichan
%========
if isfield(in,'ichan')
  if ~isempty(in.ichan)
    in.ichan = round(in.ichan); % force to integer
    if any(in.ichan(:)<1 | in.ichan(:)>size(in.tdata,2))
      tmpstr = 'Each entry in ichan(:) must be an integer ';
      tmpstr = [tmpstr,'between 1 and 2nd dim. of tdata.'];
      error(tmpstr);
    end
    in.tdata = in.tdata(:,in.ichan);
  
    if ~in.myoptions(1)
      disp('The following channels are used for modal identification:');
      disp(num2str(in.ichan(:).'));
    end
  end
  % remove ichan from IN since tdata is already modified
  in = rmfield(in,'ichan');
end

% iband
%=======
if isfield(in,'iband')
  if ~isempty(in.iband)
    in.iband = round(in.iband); % force to integer
    if any(in.iband(:)<1 | in.iband(:)>size(in.f1f2,1))
      tmpstr = 'Each entry in iband(:) must be an integer ';
      tmpstr = [tmpstr,'between 1 and 1st dim. of f1f2.'];
      error(tmpstr);
    end
    in.f0 = in.f0(in.iband);
    in.f1f2 = in.f1f2(in.iband,:);
  end
  % remove ichan from IN since tdata is already modified
  in = rmfield(in,'iband');
end

% fs
%=====
if length(in.fs(:))>1
  warning('Only fs(1) will be used.');
end
in.fs = in.fs(1);

% f1f2
%======
if ~isnumeric(in.f1f2)
  error('f1f2 should be a numerical array.');
end
if size(in.f1f2,2)~=2
  error('f1f2 should have 2 columns')
end
nb = size(in.f1f2,1); % no. of freq. bands

% f0
%======
if ~iscell(in.f0)
  error('f0 should be a cell.');
end
if ~isequal(size(in.f0),[nb 1])
  error('The 1st dim. of f0 and f1f2 should be equal.');
end

% nseg
%======
if rem(in.nseg,1) || (in.nseg<1) || (in.nseg>size(in.tdata,1))
  error('nseg must be an integer betwen 1 and first dim. of tdata');
end

% noleak does not work for nseg>1
%=================================
if in.nseg>1 && in.noleak
  error('Bias correction algorithm (NOLEAK==1) does not work for NSEG>1');
end

%==========================================================================  
% run through different bands
%==========================================================================  

o = cell(1,nb);
in0 = in;
for ii = 1:nb
  in = in0;
  in.f0 = in0.f0(ii);
  in.f1f2 = in0.f1f2(ii,:);
  o{ii} = feval(fun,in);
end

% combine fields in different bands
out = o{1};
varx = fieldnames(o{1}).'; % cell array of field names
for ii=2:nb
  for xx = varx
    if isstruct(out.(xx{:}))
      vary = fieldnames(out.(xx{:})).';
      for yy=vary
        out.(xx{:}).(yy{:}) = [out.(xx{:}).(yy{:}),o{ii}.(xx{:}).(yy{:})];
      end
    else
      out.(xx{:}) = [out.(xx{:}),o{ii}.(xx{:})];
    end
  end
end

% arrange field orders
varx = {'f','z','Sii','Se','phi','rms','sn','k','coefv','coefv0','coefv1',...
        'S','time_req','cmat','loglik'};
% varx = {'f','z','Sii','Se','phi','rms','sn','k','coefv','coefv0','coefv1',...
%         'S','time_req','cmat','ev','d','ff'};
out = orderfields(out,varx);

