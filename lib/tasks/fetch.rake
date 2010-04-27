namespace :fetch do

  desc "Fetch SCOP parseable files from MRC-LMB Web site"
  task :scop => [:environment] do

    require "open-uri"
    require "hpricot"

    scop_uri  = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/"
    mink_scop_dir   = configatron.mink_scop_dir
    links     = Hash.new(0)

    Hpricot(open(scop_uri)).search("//a") do |link|
      if link['href'] && link['href'] =~ /(dir\S+)\_(\S+)/
        stem, version = $1, $2.to_f
        links[stem] = version if links[stem] < version
      end
    end

    links.each do |stem, version|
      link = "#{stem}_#{version}"
      File.open(mink_scop_dir.join(link), 'w') do |f|
        f.puts open(scop_uri + "/#{link}").read
        $logger.info "Downloading #{link}: done"
      end
    end
  end


  desc "Fetch NR95 SCOP PDB-style files from local mirror"
  task :scop_pdb => [:environment] do

    local_scop_dir      = configatron.scop_dir
    local_scop_pdb_dir  = local_scop_dir.join("pdbstyle")
    local_scop_seq_dir  = local_scop_dir.join("scopseq")
    local_astral95      = local_scop_seq_dir.join("astral-scopdom-seqres-gd-sel-gs-bib-95-1.75.fa")
    mink_scop_pdb_dir   = configatron.mink_scop_pdb_dir

    mkdir_p mink_scop_pdb_dir

    IO.foreach(local_astral95) do |line|
      if line =~ %r{^>(\S+)\s+}
        sid = $1
        dom_sid = sid.gsub(%r{^g}, 'd')
        dom_pdb = local_scop_pdb_dir.join(dom_sid[2..3], "#{dom_sid}.ent")

        if !File.exists? dom_pdb
          $logger.error "Cannot find #{dom_pdb}"
          exit 1
        end
        cp dom_pdb, mink_scop_pdb_dir.join("#{dom_sid}.pdb")
      end
    end
  end

end
