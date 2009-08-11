class MinkSearchHit < ActiveRecord::Base

  belongs_to :mink_search

  belongs_to :norm_mink_vector

end
