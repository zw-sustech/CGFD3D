all:
	mex -g -v -f ./mexopts.R2008a.sh -output vanilla_mexnc mexgateway.c netcdf2.c netcdf3.c common.c 
	mv *.mex*  ../private

clean:
	rm -f *.o mexnc.mex*
