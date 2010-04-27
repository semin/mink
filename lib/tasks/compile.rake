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

end
