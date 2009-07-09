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
          logger.info ">>> Downloading #{link}: done"
        end
      end
    end

  end
end
