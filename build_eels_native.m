function info = build_eels_native()
%BUILD_EELS_NATIVE Compile the C++ Stage1/Stage2 kernels for this MATLAB host.
root = fileparts(mfilename('fullpath'));
src = fullfile(root,'native','eels_native_mex.cpp');
clear eels_native_v22_mex
args = {'-R2018a','-O',src,'-outdir',root,'-output','eels_native_v22_mex','-lmwblas','-lmwlapack'};
if ispc
    args = [args, {'COMPFLAGS=$COMPFLAGS /std:c++17 /O2'}];
else
    args = [args, {'CXXFLAGS=$CXXFLAGS -std=c++17 -O3'}];
end
mex(args{:});
rehash;
info = struct('version',eels_native_v22_mex('version'),'matlab',version, ...
    'binary',fullfile(root,['eels_native_v22_mex.' mexext]));
fprintf('Built %s\n%s\n',info.version,info.binary);
end
