function test_varputg ( ncfile )
% TEST_VARPUTG
%
% This routine tests VARGETG, VARPUTG
%
% Test 1:  test VARPUTG/VARGETG with double precision data
% Test 2:  test VARPUTG/VARGETG with float data, should not be accepted
% Test 3:  
% Test 4:  test 1D
% Test 005:  test 2D VARPUTG/VARGETG 
% Test 006:  test 2D VARPUTG/VARGETG on lower right quadrant 
% Test 010:  test writing a short datum to a double precision variable
% Test 100:  VARPUTG with a bad ncid
% Test 101:  VARGETG with a bad ncid
% Test 102:  VARPUTG with a bad varid
% test_minus_one_count:  count argument is negative

if ( nargin < 1 )
	ncfile = 'foo.nc';
end

mexnc ( 'setopts', 0 );

create_testfile ( ncfile );
test_double_precision ( ncfile );
test_read_single      ( ncfile );
test_read_short       ( ncfile );

test_neg_float_input ( ncfile );

test_003 ( ncfile );
test_004 ( ncfile );
test_005 ( ncfile );
test_006 ( ncfile );
test_010 ( ncfile );
test_100 ( ncfile );
test_101 ( ncfile );
test_102 ( ncfile );


fprintf ( 1, 'VARPUTG succeeded\n' );
fprintf ( 1, 'VARGETG succeeded\n' );

return




%--------------------------------------------------------------------------
function create_testfile ( ncfile )


%
% ok, first create this baby.
[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end



%
% Create the fixed dimension.  
len_x = 100;
len_y = 200;
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', len_x );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[ydimid, status] = mexnc ( 'def_dim', ncid, 'y', len_y );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[z_double_varid, status] = mexnc ( 'def_var', ncid, 'z_double', nc_double, 1, xdimid ); %#ok<ASGLU>
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[z_float_varid, status] = mexnc ( 'def_var', ncid, 'z_float', nc_float, 1, xdimid ); %#ok<ASGLU>
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end



[z_short_varid, status] = mexnc ( 'def_var', ncid, 'z_short', nc_short, 1, xdimid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[twod_varid, status] = mexnc ( 'def_var', ncid, 'twoD', nc_double, 2, [ydimid xdimid] ); %#ok<ASGLU>
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

eps = 0.01;
status = mexnc ( 'put_att_double', ncid, z_short_varid, 'scale_factor', nc_double, 1, eps );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


status = mexnc ( 'put_att_double', ncid, z_short_varid, 'add_offset', nc_double, 1, 0.00 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[status] = mexnc ( 'enddef', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return















%--------------------------------------------------------------------------
function test_read_short ( ncfile )
% a short dataset should return double precision data


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'z_short');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = 1:50;

status = mexnc ( 'VARPUTG', ncid, varid, 0, 25, 2, [], input_data(1:2:end) );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

status = mexnc ( 'VARPUTG', ncid, varid, 1, 25, 2, [], input_data(2:2:end) );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[output_data1, status] = mexnc ( 'VARGETG', ncid, varid,0,50,1);
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

if ~strcmp(class(output_data1),'double')
	error('data read by vargetg should be double precision');
end

d = max(abs(output_data1(:)-input_data(:)))';
if (any(d))
	error ( 'values written by VARGETG were not as expected.\n'  );
end






% Now get strided output
[output_data2, status] = mexnc ( 'VARGETG', ncid, varid, 0, 25, 2 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data2 = output_data2(:);

d = max(abs(output_data2-double(single(input_data(1:2:50)'))))';
if (any(d))
	error ( 'values written by VARGETG were not as expected\n'  );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return






%--------------------------------------------------------------------------
function test_read_single ( ncfile )
% a single precision dataset should return double precision data


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'z_float');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = 1:50;

status = mexnc ( 'VARPUTG', ncid, varid, 0, 25, 2, [], input_data(1:2:end) );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

status = mexnc ( 'VARPUTG', ncid, varid, 1, 25, 2, [], input_data(2:2:end) );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[output_data1, status] = mexnc ( 'VARGETG', ncid, varid, 0, 50, 1 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

if ~strcmp(class(output_data1),'double')
	error('data read by vargetg should be double precision');
end
output_data1 = output_data1(:);

d = max(abs(output_data1-double(single(input_data(:)))))';
if (any(d))
	error ( 'values written by VARGETG were not as expected.\n'  );
end






% Now get strided output
[output_data2, status] = mexnc ( 'VARGETG', ncid, varid, 0, 25, 2 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data2 = output_data2(:);

d = max(abs(output_data2-double(single(input_data(1:2:50)'))))';
if (any(d))
	error ( 'values written by VARGETG were not as expected\n'  );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return







%--------------------------------------------------------------------------
function test_double_precision ( ncfile )


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = 1:50;

status = mexnc ( 'VARPUTG', ncid, z_double_varid, 0, 25, 2, [], input_data(1:2:end) );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

status = mexnc ( 'VARPUTG', ncid, z_double_varid, 1, 25, 2, [], input_data(2:2:end) );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[output_data1, status] = mexnc ( 'VARGETG', ncid, z_double_varid, 0, 50, 1 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data1 = output_data1(:);

d = max(abs(output_data1-input_data(:)))';
if (any(d))
	error ( 'values written by VARGET do not match what was retrieved by VARPUT\n'  );
end






% Now get strided output
[output_data2, status] = mexnc ( 'VARGETG', ncid, z_double_varid, 0, 25, 2 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data2 = output_data2(:);

d = max(abs(output_data2-input_data(1:2:50)'))';
if (any(d))
	error ( 'values written by VARGET do not match what was retrieved by VARPUT\n'  );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return








%--------------------------------------------------------------------------
function test_neg_float_input ( ncfile )
% varputg should only accept double or char input


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'z_float');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = single(rand(25,1));

try
	mexnc ( 'VARPUTG', ncid, varid, 1, 25, 2, [], input_data );
	fail = true;
catch %#ok<CTCH>
    fail = false;
end;
if fail
	error ( 'Succeeded when it should have failed' );
end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return









%--------------------------------------------------------------------------
function test_010 ( ncfile )


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = int16(rand(25,1)*100);
try
	mexnc ( 'VARPUTG', ncid, varid, 1, 25, 2, [], input_data );
    fail = true;
catch %#ok<CTCH>
    fail = false;
end
if fail
	error ( 'succeeded when it should have failed' );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return









%--------------------------------------------------------------------------
function test_100 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(25,1);
status = mexnc ( 'VARPUTG', -100, z_double_varid, 1, 25, 2, [], input_data );
if status ~= -1
	error( 'VARPUTG succeeded with a bad ncid' );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return







%--------------------------------------------------------------------------
function test_101 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[data,status] = mexnc ( 'VARGETG', -100, z_double_varid, 1, 25, 2 ); %#ok<ASGLU>
if status ~= -1
	error('VARGET succeeded with a bad ncid');
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return








%--------------------------------------------------------------------------
function test_102 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(25,1);
status = mexnc ( 'VARPUTG', ncid, -500, 1, 25, 2, [], input_data );
if status ~= -1
	error('VARPUTG succeeded with a bad varid');
end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return











%--------------------------------------------------------------------------
function test_003 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_short_varid, status] = mexnc('INQ_VARID', ncid, 'z_short');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[scale_factor, status] = mexnc('GET_ATT_DOUBLE', ncid, z_short_varid, 'scale_factor');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(100,1);
input_data = input_data(1:2:end);
r = size(input_data,1);
status = mexnc ( 'VARPUTG', ncid, z_short_varid, 0, r, 2, [], input_data', 1 );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	error( '%s:  VARPUTG failed, (%s)\n', ncerr_msg );
end


[output_data, status] = mexnc ( 'VARGETG', ncid, z_short_varid, 0, r, 2, [], 1 );
if ( status ~= 0 )
	error ( '%s:  ''%s''\n', mfilename,  mexnc ( 'strerror', status ) );
end

output_data = output_data(:);

d = max(abs(output_data-input_data))';
ind = find ( d > scale_factor/2 );
if (any(ind))
	error( 'values written by VARPUTG do not match what was retrieved by VARGETG\n'  );
end



status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


%--------------------------------------------------------------------------
function test_004 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_short_varid, status] = mexnc('INQ_VARID', ncid, 'z_short');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end



input_data = rand(100,1)*10;
input_data = input_data(1:2:end);
input_data = floor(input_data);
r = size(input_data,1);
status = mexnc ( 'VARPUTG', ncid, z_short_varid, 0, r, 2, [], input_data', 0 );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	error( '%s:  VARPUTG failed, (%s)\n', mfilename, ncerr_msg );
end


[output_data, status] = mexnc ( 'VARGETG', ncid, z_short_varid, 0, r, 2, [], 0 );
if ( status ~= 0 )
	error ( '%s:  ''%s''\n', mfilename,  mexnc ( 'strerror', status ) );
end

output_data = output_data(:);

d = max(abs(output_data-input_data))';
ind = find ( d > 0 );
if (any(ind))
	error ( 'values written by VARPUTG do not match what was retrieved by VARGETG\n'  );
end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


%--------------------------------------------------------------------------
function test_005 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'twoD');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[xdimid, status] = mexnc('INQ_DIMID', ncid, 'x');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[ydimid, status] = mexnc('INQ_DIMID', ncid, 'y');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[len_x, status] = mexnc('INQ_DIMLEN', ncid, xdimid);
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[len_y, status] = mexnc('INQ_DIMLEN', ncid, ydimid);
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = 1:len_y*len_x;
input_data = reshape(input_data,[len_y len_x]);
input_data = input_data(1:2:end,1:3:end);

status = mexnc ( 'VARPUTG', ncid, varid, [0 0], size(input_data), [2 3], [], input_data' );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	error ( '%s:  VARPUT failed, (%s)\n', mfilename, ncerr_msg );
end


[output_data, status] = mexnc ( 'VARGETG', ncid, varid, [0 0], size(input_data), [2 3] );
if ( status ~= 0 )
	error ( '%s:  VARGET failed, msg ''%s''\n', mfilename, mexnc ( 'strerror', status ) );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data = output_data';

if (~strcmp(class(output_data),'double'))
	error ( 'data was not double precision' );
end

d = max(abs(output_data-input_data))';
ind = find ( d > 0 );
if (any(ind))
	error ( 'values written by VARPUT do not match what was retrieved by VARGET' );
end



%--------------------------------------------------------------------------
function test_006 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'twoD');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end




input_data = [2 3; 4 5; 6 7; 8 9];

status = mexnc ( 'VARPUTG', ncid, varid, [190 97], [4 2], [3 2], [], input_data' );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	error ( '%s:  VARPUT failed, (%s)\n', mfilename, ncerr_msg );
end


[output_data, status] = mexnc ( 'VARGETG', ncid, varid, [190 97], [4 2], [3 2], []);
if ( status ~= 0 )
	error ( '%s:  VARGET failed, msg ''%s''\n', mfilename, mexnc ( 'strerror', status ) );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data = output_data';

if (~strcmp(class(output_data),'double'))
	error( 'data was not double precision' );
end

d = max(abs(output_data-input_data))';
ind = find ( d > 0 );
if (any(ind))
	error ( 'values written by VARPUT do not match what was retrieved by VARGET' );
end



