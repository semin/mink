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

      vec_file = configatron.mink_dir.join("results", "vectors.dat")

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
          $logger.warn "!!! Cannot recognize this line: #{line.chomp}"
          next
        end
      end
      $logger.info ">>> Importing Minkowski vectors: done"
    end


    desc "Import Normalized Minkowski vectors"
    task :norm_mink_vectors => [:environment] do

      # find minimum values for each column of mink_vectors table
      min_area_a    = MinkVector.minimum(:area_a)
      min_r_half_a  = MinkVector.minimum(:r_half_a)
      min_std_a     = MinkVector.minimum(:std_a)
      min_area_p    = MinkVector.minimum(:area_p)
      min_r_half_p  = MinkVector.minimum(:r_half_p)
      min_std_p     = MinkVector.minimum(:std_p)
      min_mean      = MinkVector.minimum(:mean)
      min_std_mb    = MinkVector.minimum(:std_mb)
      min_kurtosis  = MinkVector.minimum(:kurtosis)
      min_skewness  = MinkVector.minimum(:skewness)
      min_area_e    = MinkVector.minimum(:area_e)
      min_std_e     = MinkVector.minimum(:std_e)
      min_is        = MinkVector.minimum(:is)

      sub_mink_vectors  = []
      max_area_a        = 0.0
      max_r_half_a      = 0.0
      max_std_a         = 0.0
      max_area_p        = 0.0
      max_r_half_p      = 0.0
      max_std_p         = 0.0
      max_mean          = 0.0
      max_std_mb        = 0.0
      max_kurtosis      = 0.0
      max_skewness      = 0.0
      max_area_e        = 0.0
      max_std_e         = 0.0
      max_is            = 0.0

      MinkVector.find_each do |mink_vector|
        sub_area_a   = mink_vector.area_a   - min_area_a
        sub_r_half_a = mink_vector.r_half_a - min_r_half_a
        sub_std_a    = mink_vector.std_a    - min_std_a
        sub_area_p   = mink_vector.area_p   - min_area_p
        sub_r_half_p = mink_vector.r_half_p - min_r_half_p
        sub_std_p    = mink_vector.std_p    - min_std_p
        sub_mean     = mink_vector.mean     - min_mean
        sub_std_mb   = mink_vector.std_mb   - min_std_mb
        sub_kurtosis = mink_vector.kurtosis - min_kurtosis
        sub_skewness = mink_vector.skewness - min_skewness
        sub_area_e   = mink_vector.area_e   - min_area_e
        sub_std_e    = mink_vector.std_e    - min_std_e
        sub_is       = mink_vector.is       - min_is

        max_area_a    = (max_area_a   < sub_area_a   ? sub_area_a   : max_area_a  )
        max_r_half_a  = (max_r_half_a < sub_r_half_a ? sub_r_half_a : max_r_half_a)
        max_std_a     = (max_std_a    < sub_std_a    ? sub_std_a    : max_std_a   )
        max_area_p    = (max_area_p   < sub_area_p   ? sub_area_p   : max_area_p  )
        max_r_half_p  = (max_r_half_p < sub_r_half_p ? sub_r_half_p : max_r_half_p)
        max_std_p     = (max_std_p    < sub_std_p    ? sub_std_p    : max_std_p   )
        max_mean      = (max_mean     < sub_mean     ? sub_mean     : max_mean    )
        max_std_mb    = (max_std_mb   < sub_std_mb   ? sub_std_mb   : max_std_mb  )
        max_kurtosis  = (max_kurtosis < sub_kurtosis ? sub_kurtosis : max_kurtosis)
        max_skewness  = (max_skewness < sub_skewness ? sub_skewness : max_skewness)
        max_area_e    = (max_area_e   < sub_area_e   ? sub_area_e   : max_area_e  )
        max_std_e     = (max_std_e    < sub_std_e    ? sub_std_e    : max_std_e   )
        max_is        = (max_is       < sub_is       ? sub_is       : max_is      )

        sub_mink_vectors << NVector[
          mink_vector.id,
          mink_vector.scop_id,
          mink_vector.sid,
          mink_vector.sunid,
          mink_vector.sccs,
          sub_area_a,
          sub_r_half_a,
          sub_std_a,
          sub_area_p,
          sub_r_half_p,
          sub_std_p,
          sub_mean,
          sub_std_mb,
          sub_kurtosis,
          sub_skewness,
          sub_area_e,
          sub_std_e,
          sub_is,
        ]
      end

      sub_mink_vectors.each do |sub_mink_vector|
        NormMinkVector.create!(
          :mink_vector_id => sub_mink_vector[0],
          :scop_id        => sub_mink_vector[1],
          :sid            => sub_mink_vector[2],
          :sunid          => sub_mink_vector[3],
          :sccs           => sub_mink_vector[4],
          :area_a         => sub_mink_vector[5]  / max_area_a,
          :r_half_a       => sub_mink_vector[6]  / max_r_half_a,
          :std_a          => sub_mink_vector[7]  / max_std_a,
          :area_p         => sub_mink_vector[8]  / max_area_p,
          :r_half_p       => sub_mink_vector[9]  / max_r_half_p,
          :std_p          => sub_mink_vector[10]  / max_std_p,
          :mean           => sub_mink_vector[11] / max_mean,
          :std_mb         => sub_mink_vector[12] / max_std_mb,
          :kurtosis       => sub_mink_vector[13] / max_kurtosis,
          :skewness       => sub_mink_vector[14] / max_skewness,
          :area_e         => sub_mink_vector[15] / max_area_e,
          :std_e          => sub_mink_vector[16] / max_std_e,
          :is             => sub_mink_vector[17] / max_is
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
        :max_area_a       => MinkVector.minimum(:area_a),
        :max_r_half_a     => MinkVector.minimum(:r_half_a),
        :max_std_a        => MinkVector.minimum(:std_a),
        :max_area_p       => MinkVector.minimum(:area_p),
        :max_r_half_p     => MinkVector.minimum(:r_half_p),
        :max_std_p        => MinkVector.minimum(:std_p),
        :max_mean         => MinkVector.minimum(:mean),
        :max_std_mb       => MinkVector.minimum(:std_mb),
        :max_kurtosis     => MinkVector.minimum(:kurtosis),
        :max_skewness     => MinkVector.minimum(:skewness),
        :max_area_e       => MinkVector.minimum(:area_e),
        :max_std_e        => MinkVector.minimum(:std_e),
        :max_is           => MinkVector.minimum(:is),
        :submax_area_a    => max_area_a,
        :submax_r_half_a  => max_r_half_a,
        :submax_std_a     => max_std_a,
        :submax_area_p    => max_area_p,
        :submax_r_half_p  => max_r_half_p,
        :submax_std_p     => max_std_p,
        :submax_mean      => max_mean,
        :submax_std_mb    => max_std_mb,
        :submax_kurtosis  => max_kurtosis,
        :submax_skewness  => max_skewness,
        :submax_area_e    => max_area_e,
        :submax_std_e     => max_std_e,
        :submax_is        => max_is
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


    desc "Import mink_vector_similarities"
    task :norm_mink_vector_similarities => [:environment] do

      norm_mink_vectors = NormMinkVector.all
      norm_mink_vectors.combination(2).each do |norm_mink_vector1, norm_mink_vector2|
        dist = norm_mink_vector1.euclidean_distance_to norm_mink_vector2
        MinkVectorSimilarity.create!(:norm_mink_vector_id         => norm_mink_vector1,
                                     :similar_norm_mink_vector_id => norm_mink_vector2,
                                     :distance                    => dist)
      end
    end

  end
end
