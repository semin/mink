class GiSearchHit < ActiveRecord::Base

  belongs_to :gi_search

  belongs_to :norm_gi_vector

end
