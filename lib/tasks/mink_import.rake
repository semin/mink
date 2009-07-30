namespace :mink do
  namespace :import do

    desc "Import SCOP hierarchies and descriptions"
    task :scop => [:environment] do

      tmp_dir = Rails.root.join("tmp")

      hie_file = Dir[tmp_dir.join('*hie*scop*').to_s][0]
      des_file = Dir[tmp_dir.join('*des*scop*').to_s][0]

      # Create a hash for description of scop entries,
      # and set a description for 'root' scop entry with sunid, '0'
      scop_des      = Hash.new
      scop_des['0'] = {
        :sunid        => '0',
        :stype        => 'root',
        :sccs         => 'root',
        :sid          => 'root',
        :description  => 'root',
      }

      # # dir.des.scop.txt
      # 46456   cl      a       -       All alpha proteins [46456]
      # 46457   cf      a.1     -       Globin-like
      # 46458   sf      a.1.1   -       Globin-like
      # 46459   fa      a.1.1.1 -       Truncated hemoglobin
      # 46460   dm      a.1.1.1 -       Protozoan/bacterial hemoglobin
      # 46461   sp      a.1.1.1 -       Ciliate (Paramecium caudatum) [TaxId: 5885]
      # 14982   px      a.1.1.1 d1dlwa_ 1dlw A:
      # 100068  px      a.1.1.1 d1uvya_ 1uvy A:
      IO.foreach(des_file) do |line|
        next if line =~ /^#/ || line.blank?
        sunid, stype, sccs, sid, description = line.chomp.split(/\t/)
        sccs  = nil if sccs =~ /unassigned/
          sid   = nil if sid  =~ /unassigned/
          scop_des[sunid] = {
          :sunid        => sunid,
          :stype        => stype,
          :sccs         => sccs,
          :sid          => sid,
          :description  => description
        }
      end

      # # dir.hie.scop.txt
      # 46460   46459   46461,46462,81667,63437,88965,116748
      # 14982   46461   -
      IO.readlines(hie_file).each_with_index do |line, i|
        next if line =~ /^#/ || line.blank?

        self_sunid, parent_sunid, children_sunids = line.chomp.split(/\t/)
        current_scop = Scop.factory_create!(scop_des[self_sunid])

        unless self_sunid.to_i == 0
          parent_scop = Scop.find_by_sunid(parent_sunid)
          current_scop.move_to_child_of(parent_scop)
        end
      end
      $logger.info ">>> Importing SCOP: done"
    end # task :scop


    desc "Import Minkowski vectors"
    task :mink_vectors => [:environment] do

      vec_file = Rails.root.join("minkscop", "results", "vectors.dat")

      unless File.exists? vec_file
        $logger.error "!!! #{vec_file} doesn not exist"
        exit 1
      end

      IO.foreach(vec_file) do |line|
        columns = line.chomp.split(/\s+/)

        if columns.size == 14
          dom = Scop.find_by_sid(columns[0])

          if dom.nil?
            $logger.error "!!! Cannot find SCOP domain, #{columns[0]}"
            exit 1
          end

          dom.create_mink_vector(:sid       => columns[0],
                                 :sunid     => dom.sunid,
                                 :area_a    => columns[1],
                                 :r_half_a  => columns[2],
                                 :std_a     => columns[3],
                                 :area_p    => columns[4],
                                 :r_half_p  => columns[5],
                                 :std_p     => columns[6],
                                 :mean      => columns[7],
                                 :std_mb    => columns[8],
                                 :kurtosis  => columns[9],
                                 :skewness  => columns[10],
                                 :area_e    => columns[11],
                                 :std_e     => columns[12],
                                 :is        => columns[13])
        else
          $logger.warn "!!! Cannot recognize this line: #{line.chomp}"
          next
        end
      end
      $logger.info ">>> Importing Minkowski vectors: done"
    end

  end
end
