classdef EELS_Checkpoint
    % Append-only result blocks; publish a small manifest after each block.
    % Legacy single-file checkpoints remain readable.
    methods (Static)
        function cp = read(file)
            s=load(file,'checkpoint');cp=s.checkpoint;
            if ~isfield(cp,'result_parts'),return;end
            cp.results=cell(1,cp.completed);offset=0;
            for k=1:numel(cp.result_parts)
                part=load(fullfile(fileparts(file),cp.result_parts{k}),'results');
                n=numel(part.results);
                if offset+n>cp.completed,error('EELS:CheckpointMismatch','Invalid checkpoint block count.');end
                cp.results(offset+(1:n))=part.results;offset=offset+n;
            end
            if offset~=cp.completed,error('EELS:CheckpointMismatch','Incomplete checkpoint blocks.');end
        end

        function state = write(settings,np,results,times,pixels,completed,state)
            if nargin<7||isempty(state),state=struct('saved',0,'parts',{{}});end
            if completed<state.saved,error('EELS:CheckpointMismatch','Checkpoint prefix moved backwards.');end
            folder=settings.test_output_dir;
            parts=fullfile(folder,'checkpoint_parts');if ~isfolder(parts),mkdir(parts);end
            if completed>state.saved
                results=results(state.saved+1:completed); %#ok<NASGU>
                partfile=[tempname(parts),'.mat'];
                save(partfile,'results','-v7.3');
                [~,name,ext]=fileparts(partfile);
                state.parts{end+1}=fullfile('checkpoint_parts',[name,ext]);
            end
            checkpoint=struct('format','eels_incremental_v1','completed',completed, ...
                'total',np,'result_parts',{state.parts},'pixel_time',times(1:completed), ...
                'pixels',pixels(1:completed), ...
                'settings',EELS_PearsonSparse_FullFit.strip_runtime_settings(settings)); %#ok<NASGU>
            temporary=[tempname(folder),'.mat'];
            save(temporary,'checkpoint','-v7');
            [ok,msg]=movefile(temporary,fullfile(folder,'datacube_fit_checkpoint.mat'),'f');
            if ~ok,error('EELS:CheckpointWrite','%s',msg);end
            state.saved=completed;
        end
    end
end
