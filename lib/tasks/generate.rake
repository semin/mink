namespace :generate do

  desc "Generate a figure for each SCOP domain only"
  task :domsolofig => [:environment] do

    cwd = pwd
    dir = configatron.scop_fig_dir
    if File.exists? dir
      rm_rf dir
    end
    mkdir_p dir
    chdir dir

    mink_vectors = MinkVector.all
    mink_vectors.each_with_index do |mink_vector, i|
      domain    = mink_vector.scop_domain
      dom_sid   = domain.sid.gsub(/^g/, "d")
      dom_sunid = domain.sunid
      dom_pdb   = configatron.scop_pdb_dir.join(dom_sid[2..3], "#{dom_sid}.ent")

      if !File.size? dom_pdb
        $logger.error "!!! Cannot find #{dom_pdb}"
        next
      end

      cp dom_pdb, "."

      pdb     = File.basename(dom_pdb)
      stem    = File.basename(pdb, ".ent")
      input   = "#{stem}.molinput"
      fig500  = "#{stem}-solo-500.png"
      fig100  = "#{stem}-solo-100.png"

      if File.size?(fig500) && File.size?(fig100)
        $logger.warn "!!! Skipped SCOP domain, #{stem}, figures are already created"
        next
      end

      mol_input       = `molauto -notitle -nice #{pdb}`.split("\n")
      mol_input[5,0]  = "  background grey 1;"

      File.open(input, "w") { |f| f.puts mol_input.join("\n") }
      system "molscript -r < #{input} | render -png #{fig500} -size500x500"
      system "convert #{fig500} -resize 100x100 #{fig100}"

      rm pdb
      rm input
    end

    chdir cwd
  end


  desc "Generate a MINK vector table file in CSV format"
  task :mink_vectors_csv => [:environment] do

    csv = configatron.mink_vectors_csv

    File.open(csv, 'w') do |file|
      MinkVector.find_each do |mink_vector|
        file.puts [
          mink_vector.id,
          mink_vector.area_a,
          mink_vector.r_half_a,
          mink_vector.std_a,
          mink_vector.area_p,
          mink_vector.r_half_p,
          mink_vector.std_p,
          mink_vector.mean,
          mink_vector.std_mb,
          mink_vector.kurtosis,
          mink_vector.skewness,
          mink_vector.area_e,
          mink_vector.std_e,
          mink_vector.is
        ].join(',')
      end
    end
    $logger.info "Generating #{csv}: done"
  end


  desc "Generate a MINK vector similarity table file in CSV format"
  task :mink_vector_similarities_csv => [:environment] do

    csv = configatron.mink_vectors_csv
    dmp = configatron.mink_vector_similarities_csv
    cmd = [
      configatron.calculate_distances_bin,
      "13", "<",
      csv.to_s,
      ">", dmp
    ].join(' ')

    sh cmd
    $logger.info "Generating #{dmp}: done"
  end


  desc "Generate a normalized MINK vector table file in CSV format"
  task :norm_mink_vectors_csv => [:environment] do

    csv = configatron.norm_mink_vectors_csv

    File.open(csv, 'w') do |file|
      NormMinkVector.find_each do |norm_mink_vector|
        file.puts [
          norm_mink_vector.id,
          norm_mink_vector.area_a,
          norm_mink_vector.r_half_a,
          norm_mink_vector.std_a,
          norm_mink_vector.area_p,
          norm_mink_vector.r_half_p,
          norm_mink_vector.std_p,
          norm_mink_vector.mean,
          norm_mink_vector.std_mb,
          norm_mink_vector.kurtosis,
          norm_mink_vector.skewness,
          norm_mink_vector.area_e,
          norm_mink_vector.std_e,
          norm_mink_vector.is
        ].join(',')
      end
    end
    $logger.info "Generating #{csv}: done"
  end


  desc "Generate a normalized MINK vector similarity table file in CSV format"
  task :norm_mink_vector_similarities_csv => [:environment] do

    csv = configatron.norm_mink_vectors_csv
    dmp = configatron.norm_mink_vector_similarities_csv
    cmd = [
      configatron.calculate_distances_bin,
      "13", "<",
      csv.to_s,
      ">", dmp
    ].join(' ')

    sh cmd
    $logger.info "Generating #{dmp}: done"
  end


  desc "Generate a GI vector table file in CSV format"
  task :gi_vectors_csv => [:environment] do

    csv = configatron.gi_vectors_csv

    File.open(csv, 'w') do |file|
      GiVector.find_each do |gi_vector|
        file.puts [
          gi_vector.id,
          gi_vector.length,
          gi_vector.int12,
          gi_vector.inta12,
          gi_vector.int12_34,
          gi_vector.inta12_34,
          gi_vector.int12_a34,
          gi_vector.inta12_a34,
          gi_vector.int13_24,
          gi_vector.inta13_24,
          gi_vector.int13_a24,
          gi_vector.inta13_a24,
          gi_vector.int14_23,
          gi_vector.inta14_23,
          gi_vector.int14_a23,
          gi_vector.inta14_a23,
          gi_vector.int12_34_56,
          gi_vector.int12_35_46,
          gi_vector.int12_36_45,
          gi_vector.int13_24_56,
          gi_vector.int13_25_46,
          gi_vector.int13_26_45,
          gi_vector.int14_23_56,
          gi_vector.int14_25_36,
          gi_vector.int14_26_35,
          gi_vector.int15_23_46,
          gi_vector.int15_24_36,
          gi_vector.int15_26_34,
          gi_vector.int16_23_45,
          gi_vector.int16_24_35,
          gi_vector.int16_25_34
        ].join(',')
      end
    end
    $logger.info "Generating #{csv}: done"
  end


  desc "Generate a GI vector similarity table file in CSV format"
  task :gi_vector_similarities_csv => [:environment] do

    csv = configatron.gi_vectors_csv
    dmp = configatron.gi_vector_similarities_csv
    cmd = [
      configatron.calculate_distances_bin,
      "30", "<",
      csv.to_s,
      ">", dmp
    ].join(' ')

    sh cmd
    $logger.info "Generating #{dmp}: done"
  end


  desc "Generate a normalized GI vector table file in CSV format"
  task :norm_gi_vectors_csv => [:environment] do

    csv = configatron.norm_gi_vectors_csv

    File.open(csv, 'w') do |file|
      NormGiVector.find_each do |norm_gi_vector|
        file.puts [
          norm_gi_vector.id,
          norm_gi_vector.length,
          norm_gi_vector.int12,
          norm_gi_vector.inta12,
          norm_gi_vector.int12_34,
          norm_gi_vector.inta12_34,
          norm_gi_vector.int12_a34,
          norm_gi_vector.inta12_a34,
          norm_gi_vector.int13_24,
          norm_gi_vector.inta13_24,
          norm_gi_vector.int13_a24,
          norm_gi_vector.inta13_a24,
          norm_gi_vector.int14_23,
          norm_gi_vector.inta14_23,
          norm_gi_vector.int14_a23,
          norm_gi_vector.inta14_a23,
          norm_gi_vector.int12_34_56,
          norm_gi_vector.int12_35_46,
          norm_gi_vector.int12_36_45,
          norm_gi_vector.int13_24_56,
          norm_gi_vector.int13_25_46,
          norm_gi_vector.int13_26_45,
          norm_gi_vector.int14_23_56,
          norm_gi_vector.int14_25_36,
          norm_gi_vector.int14_26_35,
          norm_gi_vector.int15_23_46,
          norm_gi_vector.int15_24_36,
          norm_gi_vector.int15_26_34,
          norm_gi_vector.int16_23_45,
          norm_gi_vector.int16_24_35,
          norm_gi_vector.int16_25_34
        ].join(',')
      end
    end
    $logger.info "Generating #{csv}: done"
  end


  desc "Generate a normalized GI vector similarity table file in CSV format"
  task :norm_gi_vector_similarities_csv => [:environment] do

    csv = configatron.norm_gi_vectors_csv
    dmp = configatron.norm_gi_vector_similarities_csv
    cmd = [
      configatron.calculate_distances_bin,
      "30", "<",
      csv.to_s,
      ">", dmp
    ].join(' ')

    sh cmd
    $logger.info "Generating #{dmp}: done"
  end


  desc "Generate a plot showing correlation between GI and MINK"
  task :correlation_mink_gi_csv => [:environment] do

    file = File.open("./tmp/mink-gi_distance.csv", "w")
    file.puts ["MINK", "GI"].join(',')

    NormMinkVectorSimilarity.find_each do |norm_mink_vector_similarity|
      norm_mink_vector1 = norm_mink_vector_similarity.norm_mink_vector
      norm_mink_vector2 = norm_mink_vector_similarity.norm_mink_vector_target
      domain1 = norm_mink_vector1.scop_domain
      domain2 = norm_mink_vector2.scop_domain

      if ((domain1.norm_gi_vectors.size == 1) &&
          (domain2.norm_gi_vectors.size == 1))
        gi_vector1    = domain1.norm_gi_vectors.first
        gi_vector2    = domain2.norm_gi_vectors.first
        gi_distance   = gi_vector1.euclidean_distance_to gi_vector2
        mink_distance = norm_mink_vector_similarity.distance

        file.puts [mink_distance, gi_distance].join(',')
      end
    end
    file.close
  end

end
