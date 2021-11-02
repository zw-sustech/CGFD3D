function test__create ( ncfile )
% TEST__CREATE:
%
% Tests run are open with
% Test 1:   nc_clobber_mode.  Check the initial file size.
% Test 2:   nc_noclobber_mode
% Test 3:   clobber and share and 64 bit offset
% Test 4:  share mode.  Should also clobber it.
% Test 5:  share | 64bit_offset
% Test 6:  64 bit offset.  Should also clobber it.
% Test 7:  noclobber mode.  Should not succeed.
% Test 8:  only one input, should not succeed
% Test 9:  Filename is empty
% Test 10:  mode argument not supplied
%
% Basically the tests are the same as those for CREATE except we are
% using _CREATE instead.
%
% The _CREATE routine really isn't necesary anymore in NetCDF-4.  This
% is for backwards compatibility only.

if nargin == 0
	ncfile = 'foo.nc';
end


test_no_clobber_mode(ncfile);
test_no_clobber_shared_mode(ncfile);
test_no_clobber_shared_64bit_offset_mode(ncfile);
test_shared_mode(ncfile);
test_shared_64bit_offset_mode(ncfile);
test_64bit_offset_mode(ncfile);
test_no_clobber_mode(ncfile);
test_only_one_input;
test_filename_is_empty;

fprintf ( '_CREATE succeeded.\n' );

%--------------------------------------------------------------------------
function test_no_clobber_mode(ncfile)

% There seems to be a behavior change in R2010b on windows. 
initsz = 4096;
[chunksize,ncid, status] = mexnc('_create',ncfile,nc_clobber_mode,initsz);
if status, error(mexnc('strerror',status)), end
status = mexnc ( 'close', ncid );
if status, error(mexnc('strerror',status)), end

d = dir ( ncfile );
if d.bytes ~= initsz
	error('initialsize not honored.');
end






%--------------------------------------------------------------------------
function test_no_clobber_shared_mode(ncfile)
% Test 2:   nc_noclobber_mode | nc_share_mode
mode = bitor ( nc_clobber_mode, nc_share_mode );
[chunksize, ncid, status] = mexnc ( '_create', ncfile, mode, 5000 );
if status, error(mexnc('strerror',status)), end
status = mexnc ( 'close', ncid );
if status, error(mexnc('strerror',status)), end



%--------------------------------------------------------------------------
function test_no_clobber_shared_64bit_offset_mode(ncfile)
% Test 3:   clobber and share and 64 bit offset
mode = bitor ( nc_clobber_mode, nc_share_mode );
mode = bitor ( mode, nc_64bit_offset_mode );
[chunksize, ncid, status] = mexnc ( '_create', ncfile, mode, 5000 );
if status, error(mexnc('strerror',status)), end

status = mexnc ( 'close', ncid );
if status, error(mexnc('strerror',status)), end


%--------------------------------------------------------------------------
function test_shared_mode(ncfile)
% Test 4:  share mode.  Should also clobber it.
[chunksize, ncid, status] = mexnc ( '_create', ncfile, nc_share_mode, 5000 );
if status, error(mexnc('strerror',status)), end

status = mexnc ( 'close', ncid );
if status, error(mexnc('strerror',status)), end


%--------------------------------------------------------------------------
function test_shared_64bit_offset_mode(ncfile)
% Test 5:  share | 64bit_offset

mode = bitor ( nc_share_mode, nc_64bit_offset_mode );
[chunksize, ncid, status] = mexnc ( '_create', ncfile, mode, 5000 );
if status, error(mexnc('strerror',status)), end

status = mexnc ( 'close', ncid );
if status, error(mexnc('strerror',status)), end



%--------------------------------------------------------------------------
function test_64bit_offset_mode(ncfile)
% Test 6:  64 bit offset.  Should also clobber it.
[chunksize, ncid, status] = mexnc ( '_create', ncfile, nc_64bit_offset_mode, 5000 );
if status, error(mexnc('strerror',status)), end

status = mexnc ( 'close', ncid );
if status, error(mexnc('strerror',status)), end




%--------------------------------------------------------------------------
function test_only_one_input()
% Test 8:  only one input, should not succeed.  Throws an exception, 
%          because there are way too few arguments.
return
try
	[chunksize, ncid, status] = mexnc ( '_create' );
    return
catch	
	%
end

error ( 'succeeded when it should have failed' );



%--------------------------------------------------------------------------
function test_filename_is_empty()
try
	[chunksize, ncid, status] = mexnc ( '_create', '', nc_clobber_mode, 5000 );
    return
catch
    %The
end

error('succeeded when it should have failed' );


