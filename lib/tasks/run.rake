namespace :run do

  desc "Run MINK"
  task :mink => [:environment] do
    cmd = [
      configatron.batch_mink,
      "-x #{configatron.mink}",
      "-i #{configatron.mink_scop_pdb_dir}",
      "-o #{configatron.mink_scop_mink_dir}"
    ].join(" ")

    sh cmd
  end


  desc "Run proc_mink.py"
  task :proc_mink => [:environment] do

    cwd = pwd
    cd configatron.mink_scop_mink_dir
    cmd = configatron.proc_mink
    sh cmd.to_s
    cd cwd
  end


  desc "Run GI"
  task :gi => [:environment] do

    mkdir_p configatron.mink_scop_gi_dir

    cmd = [
      configatron.gi,
      configatron.mink_scop_pdb_dir.to_s + "/", # this trailing / is critical!
      configatron.mink_scop_gi_dir.join("GI.stdout"),
      configatron.mink_scop_gi_dir.join("GI.stderr")
    ].join(" ")

    sh cmd
  end

end
