classdef EELS_TuningParallel
    % Same process pool / one BLAS thread policy as Stage1/2 full-map fitting.
    methods (Static)
        function execution=prepare(workers)
            if exist('parpool','file')~=2||~license('test','Distrib_Computing_Toolbox')
                error('EELS:TuningParallel','Parallel Computing Toolbox required; turn off parallel fitting to run serially.');
            end
            pool=EELS_Parallel.prepare(workers);
            setup=parfevalOnAll(pool,@EELS_TuningParallel.workerSetup,1,fileparts(mfilename('fullpath')));
            fetchOutputs(setup);
            execution=struct('mode','parallel_points','workers',min(workers,pool.NumWorkers), ...
                'pool_workers',pool.NumWorkers,'native_backend','cpp','blas_threads_per_worker',1, ...
                'pso_parallel',false,'pool_reused_across_candidates',true);
        end
        function v=workerSetup(root)
            addpath(root);rehash;
            s=struct('native_backend','cpp');EELS_Native.enabled(s);
            v=eels_native_v22_mex('version');
        end
        function outputs=fit(e,Y,s,T,options,label)
            np=size(Y,2);outputs=cell(1,np);
            if ~options.parallel_enabled
                for j=1:np
                    drawnow limitrate;EELS_TuningParallel.checkStop(options);
                    outputs{j}=EELS_TuningParallel.point(e,Y(:,j),s,T,j==1);
                    EELS_TuningParallel.notify(options,j,np,label);
                end
                return;
            end
            pool=gcp('nocreate');
            if isempty(pool),error('EELS:TuningParallel','Parallel pool is unavailable.');end
            futures=parallel.FevalFuture.empty;ids=[];next=1;done=0;
            limit=min([np,options.parallel_workers,pool.NumWorkers]);
            try
                while done<np
                    drawnow limitrate;EELS_TuningParallel.checkStop(options);
                    while next<=np&&numel(futures)<limit
                        futures(end+1)=parfeval(pool,@EELS_TuningParallel.point,1,e,Y(:,next),s,T,next==1); %#ok<AGROW>
                        ids(end+1)=next;next=next+1; %#ok<AGROW>
                    end
                    [finished,out]=fetchNext(futures,.1);
                    if isempty(finished),continue;end
                    outputs{ids(finished)}=out;futures(finished)=[];ids(finished)=[];
                    done=done+1;EELS_TuningParallel.notify(options,done,np,label);
                end
            catch ME
                if ~isempty(futures),cancel(futures);end
                rethrow(ME);
            end
        end
        function out=point(e,y,s,T,keepResult)
            out=struct('result',[],'valley',[],'NRMSE',NaN,'reported_T_K',NaN, ...
                'status','ok','worker_pid',feature('getpid'),'native_version','');
            try
                if EELS_Native.enabled(s),out.native_version=eels_native_v22_mex('version');end
                first=EELS_PearsonSparse_FullFit.solve_stage1(e,y,s);
                current=EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(first,s);
                out.valley=EELS_ReferenceTuner.valley(current,T);out.NRMSE=current.metrics.NRMSE;
                out.reported_T_K=current.stage2.T_fit;
                if keepResult,out.result=current;end
            catch ME
                out.status=[ME.identifier,': ',ME.message];
            end
        end
        function checkStop(options)
            if ~isempty(options.stop_callback)&&options.stop_callback()
                error('EELS:PSOStop','User requested stop.');
            end
        end
        function notify(options,done,total,label)
            if isfield(options,'point_progress_callback')&&~isempty(options.point_progress_callback)
                options.point_progress_callback(done,total,label);
            end
        end
    end
end
