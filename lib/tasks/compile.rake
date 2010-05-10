namespace :compile do

  desc "Compile ./src/calculate_distances.c"
  task :calculate_distances => [:environment] do

    cmd = [
      "gcc",
      configatron.calculate_distances_src,
      "-o",
      configatron.calculate_distances_bin,
      "-lm"
    ].join(" ")

    sh cmd
  end


  desc "Compile ./src/GI.c"
  task :gi => [:environment] do

    # gcc GI.c -lm -O3 -o GI
    cmd = [
      "gcc",
      configatron.gi_src,
      "-lm",
      "-O3",
      "-o",
      configatron.gi_bin
    ].join(" ")

    sh cmd
  end


  desc "Compile ./src/GIT.c"
  task :git => [:environment] do

    # gcc GIT.c -lm -O3 -o GIT
    cmd = [
      "gcc",
      configatron.git_src,
      "-lm",
      "-O3",
      "-o",
      configatron.git_bin
    ].join(" ")

    sh cmd
  end

end
