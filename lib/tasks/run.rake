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
      configatron.gi_bin,
      configatron.mink_scop_pdb_dir.to_s + "/", # this trailing / is critical!
      configatron.mink_scop_gi_dir.join("GI.stdout"),
      configatron.mink_scop_gi_dir.join("GI.stderr")
    ].join(" ")

    sh cmd
  end


  desc "Run GIT"
  task :git => [:environment] do

    mkdir_p configatron.mink_scop_git_dir
    pdbs  = Pathname.glob(configatron.mink_scop_pdb_dir.join("*.pdb").to_s)
    pfm   = Parallel::ForkManager.new(5)

    pdbs.each_with_index do |pdb, i|
      begin
        pfm.start(pdb) and next
        sid = pdb.basename(".pdb")
        dir = configatron.mink_scop_git_dir.join(sid)
        mkdir_p dir
        cp pdb, dir
        cmd = [
          configatron.git_bin,
          dir.to_s + "/", # this trailing / is critical!
          configatron.avg_gauss_tbl_smt_rep,
          configatron.mink_scop_git_dir.join("#{sid}.out"),
          configatron.mink_scop_git_dir.join("#{sid}.err"),
          "1>/dev/null 2>&1"
        ].join(" ")
        sh cmd
        puts "#{pdb} (#{i+1}/#{pdbs.size}): done"
      rescue
        warn "Something wrong with #{pdb}"
        pfm.finish(255)
      end
      pfm.finish
    end
    pfm.wait_all_children()

    cd configatron.mink_scop_git_dir
    sh 'ls | grep "\.out" | xargs cat > GIT.stdout'
    sh 'ls | grep "\.err" | xargs cat > GIT.stderr'
    cd Rails.root
  end

end
