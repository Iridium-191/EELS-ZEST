classdef EELS_Parallel
    % One parallel layer: independent pixels on process workers, one BLAS
    % thread per worker. Keep the pool warm across Stage1/Stage2 GUI actions.
    methods (Static)
        function n = default_workers()
            try
                n=min(8,max(1,feature('numcores')));
                cluster=parcluster('local');
                n=min(n,cluster.NumWorkers);
            catch
                n=1;
            end
        end

        function pool = prepare(workers)
            pool=gcp('nocreate');
            if ~isempty(pool) && (isa(pool,'parallel.ThreadPool') || pool.NumWorkers<workers)
                delete(pool);pool=[];
            end
            if isempty(pool),pool=parpool('local',workers);end
            persistent configured_pool
            if isempty(configured_pool)||~isvalid(configured_pool)||~isequal(configured_pool,pool)
                f=parfevalOnAll(pool,@maxNumCompThreads,0,1);
                fetchOutputs(f);
                configured_pool=pool;
            end
        end

        function [results,times] = fit_chunks(spectra,ene,pixel_chunks,chunks,ny,nx,settings,workers)
            % A bounded queue keeps Stop responsive and makes checkpoint /
            % resume effective in the GUI's parallel path as well as serial.
            pool=EELS_Parallel.prepare(workers);np=sum(cellfun(@numel,chunks));
            results=cell(1,np);times=NaN(1,np);done=0;
            file=fullfile(settings.test_output_dir,'datacube_fit_checkpoint.mat');
            if settings.datacube_resume_checkpoint && isfile(file)
                cp=EELS_Checkpoint.read(file);
                pixels=[pixel_chunks{:}];
                valid=isfield(cp,'total')&&cp.total==np&&isfield(cp,'results')&& ...
                    isfield(cp,'pixels')&&numel(cp.pixels)<=np&& ...
                    isequal(cp.pixels(:),reshape(pixels(1:numel(cp.pixels)),[],1));
                if valid && isfield(cp,'settings') && isfield(cp.settings,'datacube_input_file') && ...
                        isfield(settings,'datacube_input_file')
                    valid=strcmp(cp.settings.datacube_input_file,settings.datacube_input_file);
                end
                if ~valid,error('EELS:CheckpointMismatch','Checkpoint does not match this datacube and pixel selection.');end
                done=min(numel(cp.results),np);results(1:done)=cp.results(1:done);
                if isfield(cp,'pixel_time'),times(1:done)=cp.pixel_time(1:done);end
                fprintf('Parallel checkpoint: resuming after pixel %d/%d\n',done,np);
            end
            if done==np,return;end
            for k=1:numel(chunks)
                keep=chunks{k}>done;chunks{k}=chunks{k}(keep);
                pixel_chunks{k}=pixel_chunks{k}(keep);spectra{k}=spectra{k}(keep,:);
            end
            keep=~cellfun(@isempty,chunks);chunks=chunks(keep);pixel_chunks=pixel_chunks(keep);spectra=spectra(keep);
            worker_settings=settings;worker_settings.datacube_progress_callback=[];worker_settings.datacube_stop_callback=[];
            progress_settings=settings;progress_settings.datacube_parallel_effective_workers=workers;
            timer=tic;
            checkpoint_state=[];last_save_clock=tic;
            futures=parallel.FevalFuture.empty;task_ids=[];next=1;completed=done;last_saved=done;
            queue_limit=workers;
            if pool.NumWorkers==workers,queue_limit=2*workers;end
            try
                while true
                    drawnow limitrate;
                    if EELS_PearsonSparse_FullFit.datacube_stop_requested(settings)
                        cancel(futures);
                        prefix=find(cellfun(@isempty,results),1)-1;
                        if isempty(prefix),prefix=np;end
                        EELS_Parallel.save_prefix(settings,results,times,ny,nx,np,prefix,checkpoint_state);
                        error('EELS:DatacubePaused','Parallel fitting paused; checkpoint saved through pixel %d/%d.',prefix,np);
                    end
                    if next>numel(chunks)&&isempty(futures),break;end
                    while next<=numel(chunks)&&numel(futures)<queue_limit
                        futures(end+1)=parfeval(pool,@EELS_PearsonSparse_FullFit.fit_datacube_pixel_chunk,2, ...
                            spectra{next},ene,pixel_chunks{next},ny,nx,worker_settings,[]); %#ok<AGROW>
                        task_ids(end+1)=next;next=next+1; %#ok<AGROW>
                    end
                    [finished,res,elapsed]=fetchNext(futures,0.1);
                    if isempty(finished),continue;end
                    task=task_ids(finished);idx=chunks{task};results(idx)=res;times(idx)=elapsed;
                    completed=completed+numel(idx);futures(finished)=[];task_ids(finished)=[];
                    % Report only retrieved results: queued DataQueue messages
                    % can otherwise arrive after the final future was removed.
                    EELS_PearsonSparse_FullFit.notify_progress(progress_settings,completed,np, ...
                        toc(timer),times,'datacube parallel');
                    prefix=find(cellfun(@isempty,results),1)-1;if isempty(prefix),prefix=np;end
                    if settings.datacube_checkpoint_enabled && settings.datacube_checkpoint_every>0 && ...
                            ((prefix-last_saved>=settings.datacube_checkpoint_every && toc(last_save_clock)>=30)||completed==np)
                        checkpoint_state=EELS_Parallel.save_prefix(settings,results,times,ny,nx,np,prefix,checkpoint_state);
                        last_saved=prefix;last_save_clock=tic;
                    end
                end
            catch ME
                if ~isempty(futures),cancel(futures);end
                rethrow(ME);
            end
        end

        function state = save_prefix(settings,results,times,ny,nx,np,prefix,state)
            if nargin<8,state=[];end
            % Results carry original indices even for noncontiguous samples.
            pixels=zeros(1,np);
            for k=1:prefix
                if isfield(results{k},'pixel_index'),pixels(k)=results{k}.pixel_index;
                else,pixels(k)=sub2ind([ny,nx],results{k}.pixel_y,results{k}.pixel_x);end
            end
            previously_saved=0;if ~isempty(state),previously_saved=state.saved;end
            save_timer=tic;
            state=EELS_Checkpoint.write(settings,np,results,times,pixels,prefix,state);
            if isfield(settings,'datacube_progress_every')&&settings.datacube_progress_every>0
                fprintf('Checkpoint: %d new pixels, prefix %d/%d, %.2f s\n', ...
                    prefix-previously_saved,prefix,np,toc(save_timer));
            end
        end
    end
end
