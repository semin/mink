require 'fileutils'

class MinkSearchJob < Struct.new(:mink_search_id)

  include FileUtils

  def perform
    mink_search = MinkSearch.find(mink_search_id)
    mink_search.started_at = Time.now
    mink_search.status = 'Running MINK'
    mink_search.save!

    # generate figures
    mink_search.generate_images

    pdbpath   = Pathname.new(mink_search.pdb.path)
    stem      = pdbpath.basename(pdbpath.extname)
    minkpath  = pdbpath.dirname.join("#{stem}.vec")

    system "#{configatron.mink} #{mink_search.pdb.path} > #{minkpath}"
    minkvec = `#{configatron.minkproc} #{minkpath}`.chomp.split(/\s+/).map(&:to_f)

    mink_search.ori_area_a    = minkvec[0]
    mink_search.ori_r_half_a  = minkvec[1]
    mink_search.ori_std_a     = minkvec[2]
    mink_search.ori_area_p    = minkvec[3]
    mink_search.ori_r_half_p  = minkvec[4]
    mink_search.ori_std_p     = minkvec[5]
    mink_search.ori_mean      = minkvec[6]
    mink_search.ori_std_mb    = minkvec[7]
    mink_search.ori_kurtosis  = minkvec[8]
    mink_search.ori_skewness  = minkvec[9]
    mink_search.ori_area_e    = minkvec[10]
    mink_search.ori_std_e     = minkvec[11]
    mink_search.ori_is        = minkvec[12]

    mink_search.area_a    = (minkvec[0]  - MinkValue.first.min_area_a  ) / MinkValue.first.submax_area_a
    mink_search.r_half_a  = (minkvec[1]  - MinkValue.first.min_r_half_a) / MinkValue.first.submax_r_half_a
    mink_search.std_a     = (minkvec[2]  - MinkValue.first.min_std_a   ) / MinkValue.first.submax_std_a
    mink_search.area_p    = (minkvec[3]  - MinkValue.first.min_area_p  ) / MinkValue.first.submax_area_p
    mink_search.r_half_p  = (minkvec[4]  - MinkValue.first.min_r_half_p) / MinkValue.first.submax_r_half_p
    mink_search.std_p     = (minkvec[5]  - MinkValue.first.min_std_p   ) / MinkValue.first.submax_std_p
    mink_search.mean      = (minkvec[6]  - MinkValue.first.min_mean    ) / MinkValue.first.submax_mean
    mink_search.std_mb    = (minkvec[7]  - MinkValue.first.min_std_mb  ) / MinkValue.first.submax_std_mb
    mink_search.kurtosis  = (minkvec[8]  - MinkValue.first.min_kurtosis) / MinkValue.first.submax_kurtosis
    mink_search.skewness  = (minkvec[9]  - MinkValue.first.min_skewness) / MinkValue.first.submax_skewness
    mink_search.area_e    = (minkvec[10] - MinkValue.first.min_area_e  ) / MinkValue.first.submax_area_e
    mink_search.std_e     = (minkvec[11] - MinkValue.first.min_std_e   ) / MinkValue.first.submax_std_e
    mink_search.is        = (minkvec[12] - MinkValue.first.min_is      ) / MinkValue.first.submax_is

    hits  = []
    cnt   = 0
    total = NormMinkVector.count
    mink_search.status = 'Searching'

    NormMinkVector.find_each do |norm_mink_vector|
      dist = Math::sqrt(
        (mink_search.area_a   - norm_mink_vector.area_a  )**2 +
        (mink_search.r_half_a - norm_mink_vector.r_half_a)**2 +
        (mink_search.std_a    - norm_mink_vector.std_a   )**2 +
        (mink_search.area_p   - norm_mink_vector.area_p  )**2 +
        (mink_search.r_half_p - norm_mink_vector.r_half_p)**2 +
        (mink_search.std_p    - norm_mink_vector.std_p   )**2 +
        (mink_search.mean     - norm_mink_vector.mean    )**2 +
        (mink_search.std_mb   - norm_mink_vector.std_mb  )**2 +
        (mink_search.kurtosis - norm_mink_vector.kurtosis)**2 +
        (mink_search.skewness - norm_mink_vector.skewness)**2 +
        (mink_search.area_e   - norm_mink_vector.area_e  )**2 +
        (mink_search.std_e    - norm_mink_vector.std_e   )**2 +
        (mink_search.is       - norm_mink_vector.is      )**2
      )

      if dist <= mink_search.cutoff
        mink_search.mink_search_hits.create(:mink_search_id       => mink_search.id,
                                            :norm_mink_vector_id  => norm_mink_vector.id,
                                            :distance             => dist)
      end

      cnt += 1

      if cnt % 1000 == 0
        mink_search.progress = 100 * cnt / total
        mink_search.save!
      end
    end

    mink_search.finished_at   = Time.now
    mink_search.elapsed_time  = mink_search.finished_at - mink_search.started_at
    mink_search.status        = 'Finished'
    mink_search.save!
  end

end
