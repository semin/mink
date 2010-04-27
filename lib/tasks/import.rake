namespace :import do

  desc "Import SCOP hierarchies and descriptions"
  task :scop => [:environment] do

    mink_scop_dir = configatron.mink_scop_dir

    hie_file = Dir[mink_scop_dir.join('*hie*scop*').to_s][0]
    des_file = Dir[mink_scop_dir.join('*des*scop*').to_s][0]

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
    $logger.info "Importing SCOP: done"
  end # task :scop


  desc "Import MINK vectors"
  task :mink_vectors => [:environment] do

    vec_file = configatron.mink_scop_mink_dir.join("vectors.dat")

    unless File.exists? vec_file
      $logger.error "#{vec_file} doesn not exist"
      exit 1
    end

    IO.foreach(vec_file) do |line|
      columns = line.chomp.split(/\s+/)

      if columns.size == 14
        dom = Scop.find_by_sid(columns[0])

        if dom.nil?
          $logger.error "Cannot find SCOP domain, #{columns[0]}"
          exit 1
        end

        dom.create_mink_vector(:sid       => columns[0],
                               :sunid     => dom.sunid,
                               :sccs      => dom.sccs,
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
                               :is        => columns[13],
                               :scop_class_description        => dom.scop_class.description,
                               :scop_fold_description         => dom.scop_fold.description,
                               :scop_superfamily_description  => dom.scop_superfamily.description,
                               :scop_family_description       => dom.scop_family.description,
                               :scop_protein_description      => dom.scop_protein.description,
                               :scop_species_description      => dom.scop_species.description,
                               :scop_domain_description       => dom.description)
      else
        $logger.warn "Cannot recognize this line: #{line.chomp}"
        next
      end
    end
    $logger.info "Importing Minkowski vectors: done"
  end


  desc "Import Normalized Minkowski vectors"
  task :norm_mink_vectors => [:environment] do

    min_area_a  , max_area_a   = MinkVector.minimum(:area_a)  , MinkVector.maximum(:area_a)
    min_r_half_a, max_r_half_a = MinkVector.minimum(:r_half_a), MinkVector.maximum(:r_half_a)
    min_std_a   , max_std_a    = MinkVector.minimum(:std_a)   , MinkVector.maximum(:std_a)
    min_area_p  , max_area_p   = MinkVector.minimum(:area_p)  , MinkVector.maximum(:area_p)
    min_r_half_p, max_r_half_p = MinkVector.minimum(:r_half_p), MinkVector.maximum(:r_half_p)
    min_std_p   , max_std_p    = MinkVector.minimum(:std_p)   , MinkVector.maximum(:std_p)
    min_mean    , max_mean     = MinkVector.minimum(:mean)    , MinkVector.maximum(:mean)
    min_std_mb  , max_std_mb   = MinkVector.minimum(:std_mb)  , MinkVector.maximum(:std_mb)
    min_kurtosis, max_kurtosis = MinkVector.minimum(:kurtosis), MinkVector.maximum(:kurtosis)
    min_skewness, max_skewness = MinkVector.minimum(:skewness), MinkVector.maximum(:skewness)
    min_area_e  , max_area_e   = MinkVector.minimum(:area_e)  , MinkVector.maximum(:area_e)
    min_std_e   , max_std_e    = MinkVector.minimum(:std_e)   , MinkVector.maximum(:std_e)
    min_is      , max_is       = MinkVector.minimum(:is)      , MinkVector.maximum(:is)

    submax_area_a   = max_area_a   - min_area_a
    submax_r_half_a = max_r_half_a - min_r_half_a
    submax_std_a    = max_std_a    - min_std_a
    submax_area_p   = max_area_p   - min_area_p
    submax_r_half_p = max_r_half_p - min_r_half_p
    submax_std_p    = max_std_p    - min_std_p
    submax_mean     = max_mean     - min_mean
    submax_std_mb   = max_std_mb   - min_std_mb
    submax_kurtosis = max_kurtosis - min_kurtosis
    submax_skewness = max_skewness - min_skewness
    submax_area_e   = max_area_e   - min_area_e
    submax_std_e    = max_std_e    - min_std_e
    submax_is       = max_is       - min_is

    MinkVector.find_each do |mink_vector|
      dom = mink_vector.scop_domain

      NormMinkVector.create!(
        :mink_vector_id => mink_vector.id,
        :scop_id        => mink_vector.scop_id,
        :sid            => mink_vector.sid,
        :sunid          => mink_vector.sunid,
        :sccs           => mink_vector.sccs,
        :area_a         => (mink_vector.area_a   - min_area_a  ) / submax_area_a,
        :r_half_a       => (mink_vector.r_half_a - min_r_half_a) / submax_r_half_a,
        :std_a          => (mink_vector.std_a    - min_std_a   ) / submax_std_a,
        :area_p         => (mink_vector.area_p   - min_area_p  ) / submax_area_p,
        :r_half_p       => (mink_vector.r_half_p - min_r_half_p) / submax_r_half_p,
        :std_p          => (mink_vector.std_p    - min_std_p   ) / submax_std_p,
        :mean           => (mink_vector.mean     - min_mean    ) / submax_mean,
        :std_mb         => (mink_vector.std_mb   - min_std_mb  ) / submax_std_mb,
        :kurtosis       => (mink_vector.kurtosis - min_kurtosis) / submax_kurtosis,
        :skewness       => (mink_vector.skewness - min_skewness) / submax_skewness,
        :area_e         => (mink_vector.area_e   - min_area_e  ) / submax_area_e,
        :std_e          => (mink_vector.std_e    - min_std_e   ) / submax_std_e,
        :is             => (mink_vector.is       - min_is      ) / submax_is,
        :scop_class_description        => dom.scop_class.description,
        :scop_fold_description         => dom.scop_fold.description,
        :scop_superfamily_description  => dom.scop_superfamily.description,
        :scop_family_description       => dom.scop_family.description,
        :scop_protein_description      => dom.scop_protein.description,
        :scop_species_description      => dom.scop_species.description,
        :scop_domain_description       => dom.description
      )
    end

    MinkValue.create!(
      :min_area_a       => min_area_a,
      :min_r_half_a     => min_r_half_a,
      :min_std_a        => min_std_a,
      :min_area_p       => min_area_p,
      :min_r_half_p     => min_r_half_p,
      :min_std_p        => min_std_p,
      :min_mean         => min_mean,
      :min_std_mb       => min_std_mb,
      :min_kurtosis     => min_kurtosis,
      :min_skewness     => min_skewness,
      :min_area_e       => min_area_e,
      :min_std_e        => min_std_e,
      :min_is           => min_is,
      :max_area_a       => max_area_a,
      :max_r_half_a     => max_r_half_a,
      :max_std_a        => max_std_a,
      :max_area_p       => max_area_p,
      :max_r_half_p     => max_r_half_p,
      :max_std_p        => max_std_p,
      :max_mean         => max_mean,
      :max_std_mb       => max_std_mb,
      :max_kurtosis     => max_kurtosis,
      :max_skewness     => max_skewness,
      :max_area_e       => max_area_e,
      :max_std_e        => max_std_e,
      :max_is           => max_is,
      :submax_area_a    => submax_area_a,
      :submax_r_half_a  => submax_r_half_a,
      :submax_std_a     => submax_std_a,
      :submax_area_p    => submax_area_p,
      :submax_r_half_p  => submax_r_half_p,
      :submax_std_p     => submax_std_p,
      :submax_mean      => submax_mean,
      :submax_std_mb    => submax_std_mb,
      :submax_kurtosis  => submax_kurtosis,
      :submax_skewness  => submax_skewness,
      :submax_area_e    => submax_area_e,
      :submax_std_e     => submax_std_e,
      :submax_is        => submax_is
    )
  end


  desc "Import mink_vector_similarities"
  task :mink_vector_similarities => [:environment] do

    mink_vectors = MinkVector.all
    mink_vectors.combination(2).each do |mink_vector1, mink_vector2|
      dist = mink_vector1.euclidean_distance_to mink_vector2
      MinkVectorSimilarity.create!(:mink_vector_id          => mink_vector1,
                                   :similar_mink_vector_id  => mink_vector2,
                                   :distance                => dist)
    end
  end


  desc "Import mink_vector_similarities.csv"
  task :mink_vector_similarities_csv => [:environment] do

    csv = configatron.mink_vector_similarities_csv
    cmd = [
      "mysqlimport",
      "--host=#{Rails.configuration.database_configuration[Rails.env]['host']}",
      "--local",
      '--fields-optionally-enclosed-by=\'"\'',
      '--fields-terminated-by=,',
      '--lines-terminated-by="\n"',
      "--user=semin",
      "--password",
      Rails.configuration.database_configuration[Rails.env]["database"],
      csv
    ].join(" ")

    sh cmd
  end


  desc "Import norm_mink_vector_similarities"
  task :norm_mink_vector_similarities => [:environment] do

    norm_mink_vectors = NormMinkVector.all
    norm_mink_vectors.combination(2).each do |norm_mink_vector1, norm_mink_vector2|
      dist = norm_mink_vector1.euclidean_distance_to norm_mink_vector2
      NormMinkVectorSimilarity.create!(:norm_mink_vector_id         => norm_mink_vector1,
                                       :similar_norm_mink_vector_id => norm_mink_vector2,
                                       :distance                    => dist)
    end
  end


  desc "Import norm_mink_vector_similarities.csv"
  task :norm_mink_vector_similarities_csv => [:environment] do

    csv = configatron.norm_mink_vector_similarities_csv
    cmd = [
      "mysqlimport",
      "--host=#{Rails.configuration.database_configuration[Rails.env]['host']}",
      "--local",
      '--fields-optionally-enclosed-by=\'"\'',
      '--fields-terminated-by=,',
      '--lines-terminated-by="\n"',
      "--user=semin",
      "--password",
      Rails.configuration.database_configuration[Rails.env]["database"],
      csv
    ].join(" ")

    sh cmd
  end


  desc "Import GI vectors"
  task :gi_vectors => [:environment] do

    vec_file = configatron.mink_scop_gi_dir.join("GI.stdout")

    unless File.exists? vec_file
      $logger.error "#{vec_file} does not exist"
      exit 1
    end

    IO.foreach(vec_file) do |line|
      columns = line.chomp.split(/\s+/)

      unless columns.size == 34
        $logger.warn "Cannot recognize this line: #{line.chomp}"
        next
      end

      sid = columns[0].gsub(/^g/, 'd').gsub(/\.(pdb|ent)/, '')
      cc  = columns[1]
      dom = Scop.find_by_sid(sid)

      if dom.nil?
        $logger.error "Cannot find SCOP domain, #{sid}"
        exit 1
      end

      gv = dom.gi_vectors.find_by_sid_and_chain_code(sid, cc)

      unless gv.nil?
        $logger.warn "#{sid}-#{cc} already imported, maybe NMR structure?"
        next
      end

      dom.gi_vectors.create!(:sid                           => sid,
                             :sunid                         => dom.sunid,
                             :sccs                          => dom.sccs,
                             :chain_code                    => cc,
                             :cas                           => columns[2],
                             :cas_missing                   => columns[3],
                             :int12                         => columns[4],
                             :inta12                        => columns[5],
                             :int12_34                      => columns[6],
                             :inta12_34                     => columns[7],
                             :int12_a34                     => columns[8],
                             :inta12_a34                    => columns[9],
                             :int13_24                      => columns[10],
                             :inta13_24                     => columns[11],
                             :int13_a24                     => columns[12],
                             :inta13_a24                    => columns[13],
                             :int14_23                      => columns[14],
                             :inta14_23                     => columns[15],
                             :int14_a23                     => columns[16],
                             :inta14_a23                    => columns[17],
                             :int12_34_56                   => columns[18],
                             :int12_35_46                   => columns[19],
                             :int12_36_45                   => columns[20],
                             :int13_24_56                   => columns[21],
                             :int13_25_46                   => columns[22],
                             :int13_26_45                   => columns[23],
                             :int14_23_56                   => columns[24],
                             :int14_25_36                   => columns[25],
                             :int14_26_35                   => columns[26],
                             :int15_23_46                   => columns[27],
                             :int15_24_36                   => columns[28],
                             :int15_26_34                   => columns[29],
                             :int15_26_34                   => columns[30],
                             :int16_23_45                   => columns[31],
                             :int16_24_35                   => columns[32],
                             :int16_25_34                   => columns[33],
                             :scop_class_description        => dom.scop_class.description,
                             :scop_fold_description         => dom.scop_fold.description,
                             :scop_superfamily_description  => dom.scop_superfamily.description,
                             :scop_family_description       => dom.scop_family.description,
                             :scop_protein_description      => dom.scop_protein.description,
                             :scop_species_description      => dom.scop_species.description,
                             :scop_domain_description       => dom.description)
    end
    $logger.info "Importing GI vectors: done"
  end


  desc "Import Normalized GI vectors"
  task :norm_gi_vectors => [:environment] do

    min_int12      , max_int12       = GiVector.minimum(:int12)      , GiVector.maximum(:int12)
    min_inta12     , max_inta12      = GiVector.minimum(:inta12)     , GiVector.maximum(:inta12)
    min_int12_34   , max_int12_34    = GiVector.minimum(:int12_34)   , GiVector.maximum(:int12_34)
    min_inta12_34  , max_inta12_34   = GiVector.minimum(:inta12_34)  , GiVector.maximum(:inta12_34)
    min_int12_a34  , max_int12_a34   = GiVector.minimum(:int12_a34)  , GiVector.maximum(:int12_a34)
    min_inta12_a34 , max_inta12_a34  = GiVector.minimum(:inta12_a34) , GiVector.maximum(:inta12_a34)
    min_int13_24   , max_int13_24    = GiVector.minimum(:int13_24)   , GiVector.maximum(:int13_24)
    min_inta13_24  , max_inta13_24   = GiVector.minimum(:inta13_24)  , GiVector.maximum(:inta13_24)
    min_int13_a24  , max_int13_a24   = GiVector.minimum(:int13_a24)  , GiVector.maximum(:int13_a24)
    min_inta13_a24 , max_inta13_a24  = GiVector.minimum(:inta13_a24) , GiVector.maximum(:inta13_a24)
    min_int14_23   , max_int14_23    = GiVector.minimum(:int14_23)   , GiVector.maximum(:int14_23)
    min_inta14_23  , max_inta14_23   = GiVector.minimum(:inta14_23)  , GiVector.maximum(:inta14_23)
    min_int14_a23  , max_int14_a23   = GiVector.minimum(:int14_a23)  , GiVector.maximum(:int14_a23)
    min_inta14_a23 , max_inta14_a23  = GiVector.minimum(:inta14_a23) , GiVector.maximum(:inta14_a23)
    min_int12_34_56, max_int12_34_56 = GiVector.minimum(:int12_34_56), GiVector.maximum(:int12_34_56)
    min_int12_35_46, max_int12_35_46 = GiVector.minimum(:int12_35_46), GiVector.maximum(:int12_35_46)
    min_int12_36_45, max_int12_36_45 = GiVector.minimum(:int12_36_45), GiVector.maximum(:int12_36_45)
    min_int13_24_56, max_int13_24_56 = GiVector.minimum(:int13_24_56), GiVector.maximum(:int13_24_56)
    min_int13_25_46, max_int13_25_46 = GiVector.minimum(:int13_25_46), GiVector.maximum(:int13_25_46)
    min_int13_26_45, max_int13_26_45 = GiVector.minimum(:int13_26_45), GiVector.maximum(:int13_26_45)
    min_int14_23_56, max_int14_23_56 = GiVector.minimum(:int14_23_56), GiVector.maximum(:int14_23_56)
    min_int14_25_36, max_int14_25_36 = GiVector.minimum(:int14_25_36), GiVector.maximum(:int14_25_36)
    min_int14_26_35, max_int14_26_35 = GiVector.minimum(:int14_26_35), GiVector.maximum(:int14_26_35)
    min_int15_23_46, max_int15_23_46 = GiVector.minimum(:int15_23_46), GiVector.maximum(:int15_23_46)
    min_int15_24_36, max_int15_24_36 = GiVector.minimum(:int15_24_36), GiVector.maximum(:int15_24_36)
    min_int15_26_34, max_int15_26_34 = GiVector.minimum(:int15_26_34), GiVector.maximum(:int15_26_34)
    min_int15_26_34, max_int15_26_34 = GiVector.minimum(:int15_26_34), GiVector.maximum(:int15_26_34)
    min_int16_23_45, max_int16_23_45 = GiVector.minimum(:int16_23_45), GiVector.maximum(:int16_23_45)
    min_int16_24_35, max_int16_24_35 = GiVector.minimum(:int16_24_35), GiVector.maximum(:int16_24_35)
    min_int16_25_34, max_int16_25_34 = GiVector.minimum(:int16_25_34), GiVector.maximum(:int16_25_34)

    submax_int12       = max_int12       - min_int12
    submax_inta12      = max_inta12      - min_inta12
    submax_int12_34    = max_int12_34    - min_int12_34
    submax_inta12_34   = max_inta12_34   - min_inta12_34
    submax_int12_a34   = max_int12_a34   - min_int12_a34
    submax_inta12_a34  = max_inta12_a34  - min_inta12_a34
    submax_int13_24    = max_int13_24    - min_int13_24
    submax_inta13_24   = max_inta13_24   - min_inta13_24
    submax_int13_a24   = max_int13_a24   - min_int13_a24
    submax_inta13_a24  = max_inta13_a24  - min_inta13_a24
    submax_int14_23    = max_int14_23    - min_int14_23
    submax_inta14_23   = max_inta14_23   - min_inta14_23
    submax_int14_a23   = max_int14_a23   - min_int14_a23
    submax_inta14_a23  = max_inta14_a23  - min_inta14_a23
    submax_int12_34_56 = max_int12_34_56 - min_int12_34_56
    submax_int12_35_46 = max_int12_35_46 - min_int12_35_46
    submax_int12_36_45 = max_int12_36_45 - min_int12_36_45
    submax_int13_24_56 = max_int13_24_56 - min_int13_24_56
    submax_int13_25_46 = max_int13_25_46 - min_int13_25_46
    submax_int13_26_45 = max_int13_26_45 - min_int13_26_45
    submax_int14_23_56 = max_int14_23_56 - min_int14_23_56
    submax_int14_25_36 = max_int14_25_36 - min_int14_25_36
    submax_int14_26_35 = max_int14_26_35 - min_int14_26_35
    submax_int15_23_46 = max_int15_23_46 - min_int15_23_46
    submax_int15_24_36 = max_int15_24_36 - min_int15_24_36
    submax_int15_26_34 = max_int15_26_34 - min_int15_26_34
    submax_int15_26_34 = max_int15_26_34 - min_int15_26_34
    submax_int16_23_45 = max_int16_23_45 - min_int16_23_45
    submax_int16_24_35 = max_int16_24_35 - min_int16_24_35
    submax_int16_25_34 = max_int16_25_34 - min_int16_25_34

    GiVector.find_each do |gi_vector|
      NormGiVector.create!(
        :gi_vector_id => gi_vector.id,
        :scop_id      => gi_vector.scop_id,
        :sid          => gi_vector.sid,
        :sunid        => gi_vector.sunid,
        :sccs         => gi_vector.sccs,
        :chain_code   => gi_vector.chain_code,
        :cas          => gi_vector.cas,
        :cas_missing  => gi_vector.cas_missing,
        :int12        => (gi_vector.int12       - min_int12      ) / submax_int12,
        :inta12       => (gi_vector.inta12      - min_inta12     ) / submax_inta12,
        :int12_34     => (gi_vector.int12_34    - min_int12_34   ) / submax_int12_34,
        :inta12_34    => (gi_vector.inta12_34   - min_inta12_34  ) / submax_inta12_34,
        :int12_a34    => (gi_vector.int12_a34   - min_int12_a34  ) / submax_int12_a34,
        :inta12_a34   => (gi_vector.inta12_a34  - min_inta12_a34 ) / submax_inta12_a34,
        :int13_24     => (gi_vector.int13_24    - min_int13_24   ) / submax_int13_24,
        :inta13_24    => (gi_vector.inta13_24   - min_inta13_24  ) / submax_inta13_24,
        :int13_a24    => (gi_vector.int13_a24   - min_int13_a24  ) / submax_int13_a24,
        :inta13_a24   => (gi_vector.inta13_a24  - min_inta13_a24 ) / submax_inta13_a24,
        :int14_23     => (gi_vector.int14_23    - min_int14_23   ) / submax_int14_23,
        :inta14_23    => (gi_vector.inta14_23   - min_inta14_23  ) / submax_inta14_23,
        :int14_a23    => (gi_vector.int14_a23   - min_int14_a23  ) / submax_int14_a23,
        :inta14_a23   => (gi_vector.inta14_a23  - min_inta14_a23 ) / submax_inta14_a23,
        :int12_34_56  => (gi_vector.int12_34_56 - min_int12_34_56) / submax_int12_34_56,
        :int12_35_46  => (gi_vector.int12_35_46 - min_int12_35_46) / submax_int12_35_46,
        :int12_36_45  => (gi_vector.int12_36_45 - min_int12_36_45) / submax_int12_36_45,
        :int13_24_56  => (gi_vector.int13_24_56 - min_int13_24_56) / submax_int13_24_56,
        :int13_25_46  => (gi_vector.int13_25_46 - min_int13_25_46) / submax_int13_25_46,
        :int13_26_45  => (gi_vector.int13_26_45 - min_int13_26_45) / submax_int13_26_45,
        :int14_23_56  => (gi_vector.int14_23_56 - min_int14_23_56) / submax_int14_23_56,
        :int14_25_36  => (gi_vector.int14_25_36 - min_int14_25_36) / submax_int14_25_36,
        :int14_26_35  => (gi_vector.int14_26_35 - min_int14_26_35) / submax_int14_26_35,
        :int15_23_46  => (gi_vector.int15_23_46 - min_int15_23_46) / submax_int15_23_46,
        :int15_24_36  => (gi_vector.int15_24_36 - min_int15_24_36) / submax_int15_24_36,
        :int15_26_34  => (gi_vector.int15_26_34 - min_int15_26_34) / submax_int15_26_34,
        :int15_26_34  => (gi_vector.int15_26_34 - min_int15_26_34) / submax_int15_26_34,
        :int16_23_45  => (gi_vector.int16_23_45 - min_int16_23_45) / submax_int16_23_45,
        :int16_24_35  => (gi_vector.int16_24_35 - min_int16_24_35) / submax_int16_24_35,
        :int16_25_34  => (gi_vector.int16_25_34 - min_int16_25_34) / submax_int16_25_34,
        :scop_class_description       => gi_vector.scop_class_description,
        :scop_fold_description        => gi_vector.scop_fold_description,
        :scop_superfamily_description => gi_vector.scop_superfamily_description,
        :scop_family_description      => gi_vector.scop_family_description,
        :scop_protein_description     => gi_vector.scop_protein_description,
        :scop_species_description     => gi_vector.scop_species_description,
        :scop_domain_description      => gi_vector.scop_domain_description
      )
    end

    GiValue.create!(
      :min_int12          => min_int12,
      :min_inta12         => min_inta12,
      :min_int12_34       => min_int12_34,
      :min_inta12_34      => min_inta12_34,
      :min_int12_a34      => min_int12_a34,
      :min_inta12_a34     => min_inta12_a34,
      :min_int13_24       => min_int13_24,
      :min_inta13_24      => min_inta13_24,
      :min_int13_a24      => min_int13_a24,
      :min_inta13_a24     => min_inta13_a24,
      :min_int14_23       => min_int14_23,
      :min_inta14_23      => min_inta14_23,
      :min_int14_a23      => min_int14_a23,
      :min_inta14_a23     => min_inta14_a23,
      :min_int12_34_56    => min_int12_34_56,
      :min_int12_35_46    => min_int12_35_46,
      :min_int12_36_45    => min_int12_36_45,
      :min_int13_24_56    => min_int13_24_56,
      :min_int13_25_46    => min_int13_25_46,
      :min_int13_26_45    => min_int13_26_45,
      :min_int14_23_56    => min_int14_23_56,
      :min_int14_25_36    => min_int14_25_36,
      :min_int14_26_35    => min_int14_26_35,
      :min_int15_23_46    => min_int15_23_46,
      :min_int15_24_36    => min_int15_24_36,
      :min_int15_26_34    => min_int15_26_34,
      :min_int15_26_34    => min_int15_26_34,
      :min_int16_23_45    => min_int16_23_45,
      :min_int16_24_35    => min_int16_24_35,
      :min_int16_25_34    => min_int16_25_34,
      :max_int12          => max_int12,
      :max_inta12         => max_inta12,
      :max_int12_34       => max_int12_34,
      :max_inta12_34      => max_inta12_34,
      :max_int12_a34      => max_int12_a34,
      :max_inta12_a34     => max_inta12_a34,
      :max_int13_24       => max_int13_24,
      :max_inta13_24      => max_inta13_24,
      :max_int13_a24      => max_int13_a24,
      :max_inta13_a24     => max_inta13_a24,
      :max_int14_23       => max_int14_23,
      :max_inta14_23      => max_inta14_23,
      :max_int14_a23      => max_int14_a23,
      :max_inta14_a23     => max_inta14_a23,
      :max_int12_34_56    => max_int12_34_56,
      :max_int12_35_46    => max_int12_35_46,
      :max_int12_36_45    => max_int12_36_45,
      :max_int13_24_56    => max_int13_24_56,
      :max_int13_25_46    => max_int13_25_46,
      :max_int13_26_45    => max_int13_26_45,
      :max_int14_23_56    => max_int14_23_56,
      :max_int14_25_36    => max_int14_25_36,
      :max_int14_26_35    => max_int14_26_35,
      :max_int15_23_46    => max_int15_23_46,
      :max_int15_24_36    => max_int15_24_36,
      :max_int15_26_34    => max_int15_26_34,
      :max_int15_26_34    => max_int15_26_34,
      :max_int16_23_45    => max_int16_23_45,
      :max_int16_24_35    => max_int16_24_35,
      :max_int16_25_34    => max_int16_25_34,
      :submax_int12       => submax_int12,
      :submax_inta12      => submax_inta12,
      :submax_int12_34    => submax_int12_34,
      :submax_inta12_34   => submax_inta12_34,
      :submax_int12_a34   => submax_int12_a34,
      :submax_inta12_a34  => submax_inta12_a34,
      :submax_int13_24    => submax_int13_24,
      :submax_inta13_24   => submax_inta13_24,
      :submax_int13_a24   => submax_int13_a24,
      :submax_inta13_a24  => submax_inta13_a24,
      :submax_int14_23    => submax_int14_23,
      :submax_inta14_23   => submax_inta14_23,
      :submax_int14_a23   => submax_int14_a23,
      :submax_inta14_a23  => submax_inta14_a23,
      :submax_int12_34_56 => submax_int12_34_56,
      :submax_int12_35_46 => submax_int12_35_46,
      :submax_int12_36_45 => submax_int12_36_45,
      :submax_int13_24_56 => submax_int13_24_56,
      :submax_int13_25_46 => submax_int13_25_46,
      :submax_int13_26_45 => submax_int13_26_45,
      :submax_int14_23_56 => submax_int14_23_56,
      :submax_int14_25_36 => submax_int14_25_36,
      :submax_int14_26_35 => submax_int14_26_35,
      :submax_int15_23_46 => submax_int15_23_46,
      :submax_int15_24_36 => submax_int15_24_36,
      :submax_int15_26_34 => submax_int15_26_34,
      :submax_int15_26_34 => submax_int15_26_34,
      :submax_int16_23_45 => submax_int16_23_45,
      :submax_int16_24_35 => submax_int16_24_35,
      :submax_int16_25_34 => submax_int16_25_34
    )
  end


  desc "Import gi_vector_similarities"
  task :gi_vector_similarities => [:environment] do

    gi_vectors = GiVector.all
    gi_vectors.combination(2).each do |gi_vector1, gi_vector2|
      dist = gi_vector1.euclidean_distance_to gi_vector2
      GiVectorSimilarity.create!(:gi_vector_id          => gi_vector1,
                                 :similar_gi_vector_id  => gi_vector2,
                                 :distance              => dist)
    end
  end


  desc "Import gi_vector_similarities.csv"
  task :gi_vector_similarities_csv => [:environment] do

    csv = configatron.gi_vector_similarities_csv
    cmd = [
      "mysqlimport",
      "--host=#{Rails.configuration.database_configuration[Rails.env]['host']}",
      "--local",
      '--fields-optionally-enclosed-by=\'"\'',
      '--fields-terminated-by=,',
      '--lines-terminated-by="\n"',
      "--user=semin",
      "--password",
      Rails.configuration.database_configuration[Rails.env]["database"],
      csv
    ].join(" ")

    sh cmd
  end


  desc "Import norm_gi_vector_similarities"
  task :norm_gi_vector_similarities => [:environment] do

    norm_gi_vectors = NormGiVector.all
    norm_gi_vectors.combination(2).each do |norm_gi_vector1, norm_gi_vector2|
      dist = norm_gi_vector1.euclidean_distance_to norm_gi_vector2
      NormGiVectorSimilarity.create!(:norm_gi_vector_id         => norm_gi_vector1,
                                     :similar_norm_gi_vector_id => norm_gi_vector2,
                                     :distance                  => dist)
    end
  end


  desc "Import norm_gi_vector_similarities.csv"
  task :norm_gi_vector_similarities_csv => [:environment] do

    csv = configatron.norm_gi_vector_similarities_csv
    cmd = [
      "mysqlimport",
      "--host=#{Rails.configuration.database_configuration[Rails.env]['host']}",
      "--local",
      '--fields-optionally-enclosed-by=\'"\'',
      '--fields-terminated-by=,',
      '--lines-terminated-by="\n"',
      "--user=semin",
      "--password",
      Rails.configuration.database_configuration[Rails.env]["database"],
      csv
    ].join(" ")

    sh cmd
  end

end
