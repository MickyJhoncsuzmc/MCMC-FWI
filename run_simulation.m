function run_simulation(shot_num, shot_first_x, shot_end_x, shot_first_z, shot_end_z, f0, base_dir )

if shot_num == 1
    shot_incre_x = 0;
    shot_incre_z = 0;
else 
    shot_incre_x = (shot_end_x - shot_first_x)/(shot_num-1);
    shot_incre_z = (shot_end_z - shot_first_z)/(shot_num-1);
end 
nproc = 1;  % does not work
copy_many_shot_command = ['python generate_multi_source.py' ...
                           sprintf(' --dsx %f  --osx %f  --dsz %f --osz %f --nsx %d --f0 %f --nproc %d', ... 
                           shot_incre_x, shot_first_x, shot_incre_z, shot_first_z, shot_num, f0, nproc) ' --path ' base_dir ];
specfem2d_command = './run_this_example_PLEASE_DO_NOT_REMOVE.sh';
system('rm -r run*');
system(copy_many_shot_command);

%% submit jobs
max_concurrent_proc = min(shot_num, 20);
proc_init = zeros(max_concurrent_proc, 1);
for s_idx = 1 : shot_num
    found_slot = 0;
    while (~found_slot)
        for proc_id = 1:max_concurrent_proc
            if (proc_init(proc_id)==1)
                try
                    rc = process(proc_id).exitValue();
                    found_slot = 1;
                    break;
                catch

                end
            else
                proc_init(proc_id) = 1;
                found_slot = 1;
                break;
            end
        end
        pause(0.1);
    end
    fprintf('submitting source %d to process %d\n', s_idx, proc_id);
    cd([base_dir  sprintf('/run%04d',s_idx)]);
    runtime(proc_id) = java.lang.Runtime.getRuntime();
    process(proc_id) = runtime(proc_id).exec(specfem2d_command);
end
max_concurrent_proc = length(process);
fprintf('waiting all processes to finish...\n');
process_finished = zeros(max_concurrent_proc, 1);
while(any(process_finished==0))
    for proc_id = 1:max_concurrent_proc
        try
            rc = process(proc_id).exitValue();
            process_finished(proc_id) = 1;
        catch
            
        end
    end
    pause(0.5);
end
fprintf('****** Simulation Done ******\n');
end 


