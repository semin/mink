namespace :mink do
  namespace :fetch do

    desc "Fetch SCOP parseable files from MRC-LMB Web site"
    task :scop => [:environment] do

      require "open-uri"
      require "hpricot"

      scop_uri  = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/"
      tmp_dir   = Rails.root.join("tmp")
      links     = Hash.new(0)

      Hpricot(open(scop_uri)).search("//a") do |link|
        if link['href'] && link['href'] =~ /(dir\S+)\_(\S+)/
          stem, version = $1, $2.to_f
          links[stem] = version if links[stem] < version
        end
      end

      links.each do |stem, version|
        link = "#{stem}_#{version}"
        File.open(tmp_dir.join(link), 'w') do |f|
          f.puts open(scop_uri + "/#{link}").read
          $logger.info ">>> Downloading #{link}: done"
        end
      end
    end


    desc "Fetch NR95 SCOP PDB-style files from local mirror"
    task :scop_pdb => [:environment] do

      scop_dir = Pathname.new("/BiO/Store/SCOP")
      scop_pdb = scop_dir.join("pdbstyle")
      scop_seq = scop_dir.join("scopseq")
      astral95 = scop_seq.join("astral-scopdom-seqres-gd-sel-gs-bib-95-1.75.fa")
      mink_dir = Rails.root.join("minkscop")

      IO.foreach(astral95) do |line|
        if line =~ %r{^>(\S+)\s+}
          sid = $1
          dom_sid = sid.gsub(%r{^g}, 'd')
          dom_pdb = scop_pdb.join(dom_sid[2..3], "#{dom_sid}.ent")

          if !File.exists? dom_pdb
            $logger.error "!!! Cannot find #{dom_pdb}"
            exit 1
          end
          cp dom_pdb, mink_dir
        end
      end
    end

  end
end
