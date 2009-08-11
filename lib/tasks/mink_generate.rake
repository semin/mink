namespace :mink do
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

  end
end
