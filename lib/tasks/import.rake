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
                             :length                        => columns[4],
                             :int12                         => columns[5],
                             :inta12                        => columns[6],
                             :int12_34                      => columns[7],
                             :inta12_34                     => columns[8],
                             :int12_a34                     => columns[9],
                             :inta12_a34                    => columns[10],
                             :int13_24                      => columns[11],
                             :inta13_24                     => columns[12],
                             :int13_a24                     => columns[13],
                             :inta13_a24                    => columns[14],
                             :int14_23                      => columns[15],
                             :inta14_23                     => columns[16],
                             :int14_a23                     => columns[17],
                             :inta14_a23                    => columns[18],
                             :int12_34_56                   => columns[19],
                             :int12_35_46                   => columns[20],
                             :int12_36_45                   => columns[21],
                             :int13_24_56                   => columns[22],
                             :int13_25_46                   => columns[23],
                             :int13_26_45                   => columns[24],
                             :int14_23_56                   => columns[25],
                             :int14_25_36                   => columns[26],
                             :int14_26_35                   => columns[27],
                             :int15_23_46                   => columns[28],
                             :int15_24_36                   => columns[29],
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

    min_length     , max_length      = GiVector.minimum(:length)     , GiVector.maximum(:length)
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
    min_int16_23_45, max_int16_23_45 = GiVector.minimum(:int16_23_45), GiVector.maximum(:int16_23_45)
    min_int16_24_35, max_int16_24_35 = GiVector.minimum(:int16_24_35), GiVector.maximum(:int16_24_35)
    min_int16_25_34, max_int16_25_34 = GiVector.minimum(:int16_25_34), GiVector.maximum(:int16_25_34)

    submax_length      = max_length      - min_length
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
        :length       => (gi_vector.length      - min_length     ) / submax_length,
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
      :min_length         => min_length,
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
      :min_int16_23_45    => min_int16_23_45,
      :min_int16_24_35    => min_int16_24_35,
      :min_int16_25_34    => min_int16_25_34,
      :max_length         => max_length,
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
      :max_int16_23_45    => max_int16_23_45,
      :max_int16_24_35    => max_int16_24_35,
      :max_int16_25_34    => max_int16_25_34,
      :submax_length      => submax_length,
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


  desc "Import GIT vectors"
  task :git_vectors => [:environment] do

    vec_file = configatron.mink_scop_git_dir.join("GIT.stdout")

    unless File.exists? vec_file
      $logger.error "#{vec_file} does not exist"
      exit 1
    end

    IO.foreach(vec_file) do |line|
      columns = line.chomp.split(/\s+/)

      unless columns.size == 35
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

      gv = dom.git_vectors.find_by_sid_and_chain_code(sid, cc)

      unless gv.nil?
        $logger.warn "#{sid}-#{cc} already imported, maybe NMR structure?"
        next
      end

      dom.git_vectors.create!(:sid                          => sid,
                             :sunid                         => dom.sunid,
                             :sccs                          => dom.sccs,
                             :chain_code                    => cc,
                             :cas_missing                   => columns[2],
                             :cas                           => columns[3],
                             :cube_root_cas_19_11           => columns[4],
                             :measure1                      => columns[5],
                             :measure2                      => columns[6],
                             :measure3                      => columns[7],
                             :measure4                      => columns[8],
                             :measure5                      => columns[9],
                             :measure6                      => columns[10],
                             :measure7                      => columns[11],
                             :measure8                      => columns[12],
                             :measure9                      => columns[13],
                             :measure10                     => columns[14],
                             :measure11                     => columns[15],
                             :measure12                     => columns[16],
                             :measure13                     => columns[17],
                             :measure14                     => columns[18],
                             :measure15                     => columns[19],
                             :measure16                     => columns[20],
                             :measure17                     => columns[21],
                             :measure18                     => columns[22],
                             :measure19                     => columns[23],
                             :measure20                     => columns[24],
                             :measure21                     => columns[25],
                             :measure22                     => columns[26],
                             :measure23                     => columns[27],
                             :measure24                     => columns[28],
                             :measure25                     => columns[29],
                             :measure26                     => columns[30],
                             :measure27                     => columns[31],
                             :measure28                     => columns[32],
                             :measure29                     => columns[33],
                             :measure30                     => columns[34],
                             :scop_class_description        => dom.scop_class.description,
                             :scop_fold_description         => dom.scop_fold.description,
                             :scop_superfamily_description  => dom.scop_superfamily.description,
                             :scop_family_description       => dom.scop_family.description,
                             :scop_protein_description      => dom.scop_protein.description,
                             :scop_species_description      => dom.scop_species.description,
                             :scop_domain_description       => dom.description)
    end
    $logger.info "Importing GIT vectors: done"
  end


  desc "Import Normalized GIT vectors"
  task :norm_git_vectors => [:environment] do

    min_measure1 , max_measure1  = GitVector.minimum(:measure1 ), GitVector.maximum(:measure1 )
    min_measure2 , max_measure2  = GitVector.minimum(:measure2 ), GitVector.maximum(:measure2 )
    min_measure3 , max_measure3  = GitVector.minimum(:measure3 ), GitVector.maximum(:measure3 )
    min_measure4 , max_measure4  = GitVector.minimum(:measure4 ), GitVector.maximum(:measure4 )
    min_measure5 , max_measure5  = GitVector.minimum(:measure5 ), GitVector.maximum(:measure5 )
    min_measure6 , max_measure6  = GitVector.minimum(:measure6 ), GitVector.maximum(:measure6 )
    min_measure7 , max_measure7  = GitVector.minimum(:measure7 ), GitVector.maximum(:measure7 )
    min_measure8 , max_measure8  = GitVector.minimum(:measure8 ), GitVector.maximum(:measure8 )
    min_measure9 , max_measure9  = GitVector.minimum(:measure9 ), GitVector.maximum(:measure9 )
    min_measure10, max_measure10 = GitVector.minimum(:measure10), GitVector.maximum(:measure10)
    min_measure11, max_measure11 = GitVector.minimum(:measure11), GitVector.maximum(:measure11)
    min_measure12, max_measure12 = GitVector.minimum(:measure12), GitVector.maximum(:measure12)
    min_measure13, max_measure13 = GitVector.minimum(:measure13), GitVector.maximum(:measure13)
    min_measure14, max_measure14 = GitVector.minimum(:measure14), GitVector.maximum(:measure14)
    min_measure15, max_measure15 = GitVector.minimum(:measure15), GitVector.maximum(:measure15)
    min_measure16, max_measure16 = GitVector.minimum(:measure16), GitVector.maximum(:measure16)
    min_measure17, max_measure17 = GitVector.minimum(:measure17), GitVector.maximum(:measure17)
    min_measure18, max_measure18 = GitVector.minimum(:measure18), GitVector.maximum(:measure18)
    min_measure19, max_measure19 = GitVector.minimum(:measure19), GitVector.maximum(:measure19)
    min_measure20, max_measure20 = GitVector.minimum(:measure20), GitVector.maximum(:measure20)
    min_measure21, max_measure21 = GitVector.minimum(:measure21), GitVector.maximum(:measure21)
    min_measure22, max_measure22 = GitVector.minimum(:measure22), GitVector.maximum(:measure22)
    min_measure23, max_measure23 = GitVector.minimum(:measure23), GitVector.maximum(:measure23)
    min_measure24, max_measure24 = GitVector.minimum(:measure24), GitVector.maximum(:measure24)
    min_measure25, max_measure25 = GitVector.minimum(:measure25), GitVector.maximum(:measure25)
    min_measure26, max_measure26 = GitVector.minimum(:measure26), GitVector.maximum(:measure26)
    min_measure27, max_measure27 = GitVector.minimum(:measure27), GitVector.maximum(:measure27)
    min_measure28, max_measure28 = GitVector.minimum(:measure28), GitVector.maximum(:measure28)
    min_measure29, max_measure29 = GitVector.minimum(:measure29), GitVector.maximum(:measure29)
    min_measure30, max_measure30 = GitVector.minimum(:measure30), GitVector.maximum(:measure30)

    submax_measure1  = max_measure1  - min_measure1
    submax_measure2  = max_measure2  - min_measure2
    submax_measure3  = max_measure3  - min_measure3
    submax_measure4  = max_measure4  - min_measure4
    submax_measure5  = max_measure5  - min_measure5
    submax_measure6  = max_measure6  - min_measure6
    submax_measure7  = max_measure7  - min_measure7
    submax_measure8  = max_measure8  - min_measure8
    submax_measure9  = max_measure9  - min_measure9
    submax_measure10 = max_measure10 - min_measure10
    submax_measure11 = max_measure11 - min_measure11
    submax_measure12 = max_measure12 - min_measure12
    submax_measure13 = max_measure13 - min_measure13
    submax_measure14 = max_measure14 - min_measure14
    submax_measure15 = max_measure15 - min_measure15
    submax_measure16 = max_measure16 - min_measure16
    submax_measure17 = max_measure17 - min_measure17
    submax_measure18 = max_measure18 - min_measure18
    submax_measure19 = max_measure19 - min_measure19
    submax_measure20 = max_measure20 - min_measure20
    submax_measure21 = max_measure21 - min_measure21
    submax_measure22 = max_measure22 - min_measure22
    submax_measure23 = max_measure23 - min_measure23
    submax_measure24 = max_measure24 - min_measure24
    submax_measure25 = max_measure25 - min_measure25
    submax_measure26 = max_measure26 - min_measure26
    submax_measure27 = max_measure27 - min_measure27
    submax_measure28 = max_measure28 - min_measure28
    submax_measure29 = max_measure29 - min_measure29
    submax_measure30 = max_measure30 - min_measure30

    GitVector.find_each do |git_vector|
      NormGitVector.create!(
        :git_vector_id => git_vector.id,
        :scop_id      => git_vector.scop_id,
        :sid          => git_vector.sid,
        :sunid        => git_vector.sunid,
        :sccs         => git_vector.sccs,
        :chain_code   => git_vector.chain_code,
        :cas_missing  => git_vector.cas_missing,
        :cas          => git_vector.cas,
        :cube_root_cas_19_11 => git_vector.cube_root_cas_19_11,
        :measure1     => (git_vector.measure1  - min_measure1 ) / submax_measure1 ,
        :measure2     => (git_vector.measure2  - min_measure2 ) / submax_measure2 ,
        :measure3     => (git_vector.measure3  - min_measure3 ) / submax_measure3 ,
        :measure4     => (git_vector.measure4  - min_measure4 ) / submax_measure4 ,
        :measure5     => (git_vector.measure5  - min_measure5 ) / submax_measure5 ,
        :measure6     => (git_vector.measure6  - min_measure6 ) / submax_measure6 ,
        :measure7     => (git_vector.measure7  - min_measure7 ) / submax_measure7 ,
        :measure8     => (git_vector.measure8  - min_measure8 ) / submax_measure8 ,
        :measure9     => (git_vector.measure9  - min_measure9 ) / submax_measure9 ,
        :measure10    => (git_vector.measure10 - min_measure10) / submax_measure10,
        :measure11    => (git_vector.measure11 - min_measure11) / submax_measure11,
        :measure12    => (git_vector.measure12 - min_measure12) / submax_measure12,
        :measure13    => (git_vector.measure13 - min_measure13) / submax_measure13,
        :measure14    => (git_vector.measure14 - min_measure14) / submax_measure14,
        :measure15    => (git_vector.measure15 - min_measure15) / submax_measure15,
        :measure16    => (git_vector.measure16 - min_measure16) / submax_measure16,
        :measure17    => (git_vector.measure17 - min_measure17) / submax_measure17,
        :measure18    => (git_vector.measure18 - min_measure18) / submax_measure18,
        :measure19    => (git_vector.measure19 - min_measure19) / submax_measure19,
        :measure20    => (git_vector.measure20 - min_measure20) / submax_measure20,
        :measure21    => (git_vector.measure21 - min_measure21) / submax_measure21,
        :measure22    => (git_vector.measure22 - min_measure22) / submax_measure22,
        :measure23    => (git_vector.measure23 - min_measure23) / submax_measure23,
        :measure24    => (git_vector.measure24 - min_measure24) / submax_measure24,
        :measure25    => (git_vector.measure25 - min_measure25) / submax_measure25,
        :measure26    => (git_vector.measure26 - min_measure26) / submax_measure26,
        :measure27    => (git_vector.measure27 - min_measure27) / submax_measure27,
        :measure28    => (git_vector.measure28 - min_measure28) / submax_measure28,
        :measure29    => (git_vector.measure29 - min_measure29) / submax_measure29,
        :measure30    => (git_vector.measure30 - min_measure30) / submax_measure30,
        :scop_class_description       => git_vector.scop_class_description,
        :scop_fold_description        => git_vector.scop_fold_description,
        :scop_superfamily_description => git_vector.scop_superfamily_description,
        :scop_family_description      => git_vector.scop_family_description,
        :scop_protein_description     => git_vector.scop_protein_description,
        :scop_species_description     => git_vector.scop_species_description,
        :scop_domain_description      => git_vector.scop_domain_description
      )
    end

    GitValue.create!(
      :min_measure1     => min_measure1 ,
      :min_measure2     => min_measure2 ,
      :min_measure3     => min_measure3 ,
      :min_measure4     => min_measure4 ,
      :min_measure5     => min_measure5 ,
      :min_measure6     => min_measure6 ,
      :min_measure7     => min_measure7 ,
      :min_measure8     => min_measure8 ,
      :min_measure9     => min_measure9 ,
      :min_measure10    => min_measure10,
      :min_measure11    => min_measure11,
      :min_measure12    => min_measure12,
      :min_measure13    => min_measure13,
      :min_measure14    => min_measure14,
      :min_measure15    => min_measure15,
      :min_measure16    => min_measure16,
      :min_measure17    => min_measure17,
      :min_measure18    => min_measure18,
      :min_measure19    => min_measure19,
      :min_measure20    => min_measure20,
      :min_measure21    => min_measure21,
      :min_measure22    => min_measure22,
      :min_measure23    => min_measure23,
      :min_measure24    => min_measure24,
      :min_measure25    => min_measure25,
      :min_measure26    => min_measure26,
      :min_measure27    => min_measure27,
      :min_measure28    => min_measure28,
      :min_measure29    => min_measure29,
      :min_measure30    => min_measure30,
      :max_measure1     => max_measure1 ,
      :max_measure2     => max_measure2 ,
      :max_measure3     => max_measure3 ,
      :max_measure4     => max_measure4 ,
      :max_measure5     => max_measure5 ,
      :max_measure6     => max_measure6 ,
      :max_measure7     => max_measure7 ,
      :max_measure8     => max_measure8 ,
      :max_measure9     => max_measure9 ,
      :max_measure10    => max_measure10,
      :max_measure11    => max_measure11,
      :max_measure12    => max_measure12,
      :max_measure13    => max_measure13,
      :max_measure14    => max_measure14,
      :max_measure15    => max_measure15,
      :max_measure16    => max_measure16,
      :max_measure17    => max_measure17,
      :max_measure18    => max_measure18,
      :max_measure19    => max_measure19,
      :max_measure20    => max_measure20,
      :max_measure21    => max_measure21,
      :max_measure22    => max_measure22,
      :max_measure23    => max_measure23,
      :max_measure24    => max_measure24,
      :max_measure25    => max_measure25,
      :max_measure26    => max_measure26,
      :max_measure27    => max_measure27,
      :max_measure28    => max_measure28,
      :max_measure29    => max_measure29,
      :max_measure30    => max_measure30,
      :submax_measure1  => submax_measure1 ,
      :submax_measure2  => submax_measure2 ,
      :submax_measure3  => submax_measure3 ,
      :submax_measure4  => submax_measure4 ,
      :submax_measure5  => submax_measure5 ,
      :submax_measure6  => submax_measure6 ,
      :submax_measure7  => submax_measure7 ,
      :submax_measure8  => submax_measure8 ,
      :submax_measure9  => submax_measure9 ,
      :submax_measure10 => submax_measure10,
      :submax_measure11 => submax_measure11,
      :submax_measure12 => submax_measure12,
      :submax_measure13 => submax_measure13,
      :submax_measure14 => submax_measure14,
      :submax_measure15 => submax_measure15,
      :submax_measure16 => submax_measure16,
      :submax_measure17 => submax_measure17,
      :submax_measure18 => submax_measure18,
      :submax_measure19 => submax_measure19,
      :submax_measure20 => submax_measure20,
      :submax_measure21 => submax_measure21,
      :submax_measure22 => submax_measure22,
      :submax_measure23 => submax_measure23,
      :submax_measure24 => submax_measure24,
      :submax_measure25 => submax_measure25,
      :submax_measure26 => submax_measure26,
      :submax_measure27 => submax_measure27,
      :submax_measure28 => submax_measure28,
      :submax_measure29 => submax_measure29,
      :submax_measure30 => submax_measure30
    )
  end


  desc "Import git_vector_similarities"
  task :git_vector_similarities => [:environment] do

    git_vectors = GitVector.all
    git_vectors.combination(2).each do |git_vector1, git_vector2|
      dist = git_vector1.euclidean_distance_to git_vector2
      GitVectorSimilarity.create!(:git_vector_id          => git_vector1,
                                  :similar_git_vector_id  => git_vector2,
                                  :distance               => dist)
    end
  end


  desc "Import git_vector_similarities.csv"
  task :git_vector_similarities_csv => [:environment] do

    csv = configatron.git_vector_similarities_csv
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


  desc "Import norm_git_vector_similarities"
  task :norm_git_vector_similarities => [:environment] do

    norm_git_vectors = NormGitVector.all
    norm_git_vectors.combination(2).each do |norm_git_vector1, norm_git_vector2|
      dist = norm_git_vector1.euclidean_distance_to norm_git_vector2
      NormGitVectorSimilarity.create!(:norm_git_vector_id         => norm_git_vector1,
                                      :similar_norm_git_vector_id => norm_git_vector2,
                                      :distance                   => dist)
    end
  end


  desc "Import norm_git_vector_similarities.csv"
  task :norm_git_vector_similarities_csv => [:environment] do

    csv = configatron.norm_git_vector_similarities_csv
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
