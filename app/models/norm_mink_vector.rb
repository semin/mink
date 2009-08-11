class NormMinkVector < ActiveRecord::Base

  include Mink::Distance

  belongs_to  :mink_vector

  has_many :mink_search_hits

end
